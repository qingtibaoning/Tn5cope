#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

from tn5cope import __version__


TN5_MOSAIC_END = "CTGTCTCTTATACACATCT"
TN5_MOSAIC_END_RC = "AGATGTGTATAAGAGACAG"
TN5_MOTIFS = (
    ("TN5_ME", TN5_MOSAIC_END),
    ("TN5_ME_RC", TN5_MOSAIC_END_RC),
)
VALID_BASES = set("ACGTN")
COMP_TABLE = str.maketrans("ACGTN", "TGCAN")
ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
SEED_SIZES = (11, 9)
RESCUE_SEED_SIZES = (9, 8)
MIN_FRAGMENT_LENGTH = 20
MIN_QUERY_COVERAGE = 0.90
LOW_CONFIDENCE_IDENTITY = 0.85
TOP_HIT_MARGIN = 2
HIGH_IDENTITY_LONG = 0.95
HIGH_IDENTITY_SHORT = 0.98
RESCUE_WINDOW_SIZES = (150, 120, 90, 60, 45, 30, 24)
RESCUE_MIN_SPAN = 24
RESCUE_IDENTITY_LONG = 0.95
RESCUE_IDENTITY_SHORT = 1.0
TOP_HIT_SCORE_GAP = 2.0
JUNCTION_COORDINATE_DEFINITION = (
    "Insertion site is the 1-based reference base immediately adjacent to the "
    "Tn5-genome junction; Junction_boundary_0based is the 0-based boundary "
    "between reference bases."
)
TN5_STRUCTURE_PROFILES = (
    "unknown",
    "no_known_regulatory_elements",
    "bidirectional_outward_promoters",
    "transcriptional_terminator",
    "bidirectional_promoters_and_terminator",
)


@dataclass
class SequenceRecord:
    query_id: str
    row_number: int
    raw_sequence: str
    cleaned_sequence: str
    notes: List[str] = field(default_factory=list)


@dataclass
class MotifMatch:
    label: str
    start: int
    end: int
    mismatches: int
    observed: str


@dataclass
class FragmentCandidate:
    source: str
    sequence: str
    insertion_anchor: str
    note: str
    motif_match: Optional[MotifMatch] = None
    source_priority: int = 0


@dataclass
class Hit:
    query_id: str
    row_number: int
    candidate_source: str
    reference_seqid: str
    strand: str
    alignment_start: int
    alignment_end: int
    query_pos1_reference_coord: int
    insertion_site: int
    identity: float
    query_coverage: float
    alignment_length: int
    matched_bases: int
    comparable_bases: int
    source_priority: int = 0
    hit_rank: int = 0
    hit_count: int = 0
    hit_status: str = ""
    note: str = ""
    alignment_score: float = 0.0
    score_gap: float = 0.0
    edit_distance: int = 0
    junction_boundary_0based: int = 0
    alignment_method: str = "internal_exhaustive"


@dataclass
class GeneFeature:
    seqid: str
    start: int
    end: int
    strand: str
    gene_name: str
    gene_code: str
    locus_tag: str
    gene_id: str
    attributes: Dict[str, List[str]]


@dataclass
class CDSFeature:
    seqid: str
    start: int
    end: int
    strand: str
    gene_name: str
    locus_tag: str
    gene_id: str
    product: str
    ontology_terms: List[str]
    go_tags: List[str]
    attributes: Dict[str, List[str]]


@dataclass
class AnnotationResult:
    feature_relation: str
    gene_name: str
    gene_code: str
    putative_function: str
    ontology_term: str
    go_tags: str
    note: str
    candidate_gene_name: str = ""
    candidate_gene_code: str = ""
    candidate_relation: str = ""
    candidate_distance_bp: str = ""
    candidate_strand: str = ""
    candidate_priority: str = ""
    candidate_function: str = ""
    candidate_ontology_term: str = ""
    candidate_go_tags: str = ""
    candidate_score: str = ""
    candidate_evidence: str = ""
    candidate_operon_id: str = ""
    candidate_operon_relation: str = ""
    tn5_structure_profile: str = "unknown"


@dataclass
class NearbyGeneCandidate:
    gene: GeneFeature
    relation: str
    distance: int
    score: int
    priority: str
    evidence: List[str]
    operon_id: str = ""
    operon_relation: str = "not_assessed"


@dataclass
class CandidateEvaluation:
    candidate: FragmentCandidate
    hits: List[Hit]
    status: str
    top_hit: Optional[Hit]
    note: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="tn5cope",
        description=(
            "Locate Tn5 insertion sites, annotate gene context from GTF, "
            "and summarize genomic position clusters."
        )
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--query-csv",
        required=True,
        help="Single-column CSV containing query sequences.",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="Reference annotation in GTF format.",
    )
    parser.add_argument(
        "--genome-fasta",
        required=True,
        help="Reference genome FASTA.",
    )
    parser.add_argument(
        "--output-prefix",
        default="tn5cope_results",
        help="Prefix for output files.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Directory to place all output files in.",
    )
    parser.add_argument(
        "--aligner",
        choices=("auto", "internal", "minimap2", "blastn", "bwa"),
        default="auto",
        help=(
            "Alignment engine. auto uses minimap2, blastn, or bwa when available; "
            "otherwise it uses the built-in exhaustive ungapped fallback."
        ),
    )
    parser.add_argument(
        "--intergenic-candidate-window",
        type=int,
        default=500,
        metavar="BP",
        help=(
            "Report flanking genes within this many bp of an intergenic insertion "
            "as preliminary candidates (default: 500). This is a maximum reporting "
            "window, not a biological effect threshold. Set to 0 to disable."
        ),
    )
    parser.add_argument(
        "--operon-tsv",
        default="",
        help=(
            "Optional tab-delimited operon annotation. It must contain an operon_id "
            "column and at least one gene identifier column such as locus_tag, "
            "gene_id, gene_code, or gene_name."
        ),
    )
    parser.add_argument(
        "--tn5-structure-profile",
        choices=TN5_STRUCTURE_PROFILES,
        default="unknown",
        help=(
            "Known transcriptional structure of the inserted Tn5 construct. "
            "Unknown is conservative and makes no score adjustment."
        ),
    )
    return parser.parse_args()


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP_TABLE)[::-1]


def read_single_fasta(path: Path) -> Tuple[str, str]:
    seqid = None
    parts: List[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                new_seqid = line[1:].split()[0]
                if seqid is not None and new_seqid != seqid:
                    raise ValueError("This pipeline expects a single-contig FASTA.")
                seqid = new_seqid
                continue
            parts.append(line.upper())
    if seqid is None:
        raise ValueError(f"No FASTA record found in {path}")
    return seqid, "".join(parts)


def normalize_sequence(raw: str) -> Tuple[str, List[str]]:
    notes: List[str] = []
    seq = re.sub(r"\s+", "", raw or "").upper()
    if not seq:
        return "", notes
    invalid = sorted({base for base in seq if base not in VALID_BASES})
    if invalid:
        notes.append(f"invalid_bases_retained={''.join(invalid)}")
    ambiguous = seq.count("N")
    if ambiguous:
        notes.append(f"ambiguous_bases={ambiguous}")
    return seq, notes


def load_queries(path: Path) -> List[SequenceRecord]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        raise ValueError("Query CSV is empty.")
    first_value = re.sub(r"\s+", "", rows[0][0] if rows[0] else "").upper()
    has_header = not first_value or any(base not in VALID_BASES for base in first_value)
    data_rows = rows[1:] if has_header else rows
    records: List[SequenceRecord] = []
    for idx, row in enumerate(data_rows, start=1):
        raw = row[0] if row else ""
        cleaned, notes = normalize_sequence(raw)
        records.append(
            SequenceRecord(
                query_id=f"Q{idx:03d}",
                row_number=idx,
                raw_sequence=raw,
                cleaned_sequence=cleaned,
                notes=notes,
            )
        )
    return records


def parse_attributes(field: str) -> Dict[str, List[str]]:
    attrs: Dict[str, List[str]] = defaultdict(list)
    for key, value in ATTR_RE.findall(field):
        attrs[key].append(value)
    return dict(attrs)


def first_attr(attrs: Dict[str, List[str]], *keys: str) -> str:
    for key in keys:
        values = attrs.get(key)
        if values:
            return values[0]
    return ""


def parse_gtf(path: Path) -> Tuple[Dict[str, List[GeneFeature]], Dict[str, List[CDSFeature]]]:
    genes: Dict[str, List[GeneFeature]] = defaultdict(list)
    cdss: Dict[str, List[CDSFeature]] = defaultdict(list)
    with path.open() as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _, feature_type, start, end, _, strand, _, attrs_field = parts
            attrs = parse_attributes(attrs_field)
            start_i = int(start)
            end_i = int(end)
            if feature_type == "gene":
                gene = GeneFeature(
                    seqid=seqid,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    gene_name=first_attr(attrs, "gene"),
                    gene_code=first_attr(attrs, "old_locus_tag", "locus_tag", "gene_id"),
                    locus_tag=first_attr(attrs, "locus_tag"),
                    gene_id=first_attr(attrs, "gene_id"),
                    attributes=attrs,
                )
                genes[seqid].append(gene)
            elif feature_type == "CDS":
                go_tags = []
                for key in ("go_function", "go_process", "go_component"):
                    go_tags.extend(attrs.get(key, []))
                cds = CDSFeature(
                    seqid=seqid,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    gene_name=first_attr(attrs, "gene"),
                    locus_tag=first_attr(attrs, "locus_tag"),
                    gene_id=first_attr(attrs, "gene_id"),
                    product=first_attr(attrs, "product"),
                    ontology_terms=attrs.get("Ontology_term", []),
                    go_tags=go_tags,
                    attributes=attrs,
                )
                cdss[seqid].append(cds)
    for feature_dict in (genes, cdss):
        for seqid in feature_dict:
            feature_dict[seqid].sort(key=lambda item: (item.start, item.end))
    return dict(genes), dict(cdss)


def build_gene_lookup(genes: Dict[str, List[GeneFeature]]) -> Dict[Tuple[str, str], GeneFeature]:
    lookup: Dict[Tuple[str, str], GeneFeature] = {}
    for seqid, seq_genes in genes.items():
        for gene in seq_genes:
            for key in (gene.locus_tag, gene.gene_id, gene.gene_code, gene.gene_name):
                if key:
                    lookup[(seqid, key)] = gene
    return lookup


def normalized_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def load_operon_map(path: Optional[Path]) -> Dict[str, str]:
    """Map any supplied gene identifier to an operon identifier."""
    if path is None:
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Operon TSV has no header: {path}")
        normalized = {normalized_header(name): name for name in reader.fieldnames}
        operon_column = normalized.get("operon_id") or normalized.get("operon")
        identifier_columns = [
            normalized[key]
            for key in ("locus_tag", "gene_id", "gene_code", "gene_name", "gene")
            if key in normalized
        ]
        if not operon_column or not identifier_columns:
            raise ValueError(
                "Operon TSV must contain operon_id and at least one of: "
                "locus_tag, gene_id, gene_code, gene_name, gene."
            )
        operon_map: Dict[str, str] = {}
        for row_number, row in enumerate(reader, start=2):
            operon_id = (row.get(operon_column) or "").strip()
            if not operon_id:
                continue
            identifiers = {
                (row.get(column) or "").strip()
                for column in identifier_columns
                if (row.get(column) or "").strip()
            }
            for identifier in identifiers:
                previous = operon_map.get(identifier)
                if previous and previous != operon_id:
                    raise ValueError(
                        f"Gene identifier {identifier!r} maps to both {previous!r} "
                        f"and {operon_id!r} in {path} (line {row_number})."
                    )
                operon_map[identifier] = operon_id
    return operon_map


def operon_id_for_gene(gene: GeneFeature, operon_map: Dict[str, str]) -> str:
    for identifier in (gene.locus_tag, gene.gene_id, gene.gene_code, gene.gene_name):
        if identifier and identifier in operon_map:
            return operon_map[identifier]
    return ""


def build_kmer_index(genome: str, k: int) -> Dict[str, List[int]]:
    index: Dict[str, List[int]] = defaultdict(list)
    for pos in range(len(genome) - k + 1):
        kmer = genome[pos : pos + k]
        if "N" in kmer:
            continue
        index[kmer].append(pos)
    return dict(index)


def motif_matches(seq: str, max_mismatches: int = 2) -> List[MotifMatch]:
    matches: List[MotifMatch] = []
    for label, motif in TN5_MOTIFS:
        m_len = len(motif)
        if len(seq) < m_len:
            continue
        for start in range(len(seq) - m_len + 1):
            observed = seq[start : start + m_len]
            mismatches = sum(1 for a, b in zip(observed, motif) if a != b)
            if mismatches <= max_mismatches:
                matches.append(
                    MotifMatch(
                        label=label,
                        start=start,
                        end=start + m_len,
                        mismatches=mismatches,
                        observed=observed,
                    )
                )
    matches.sort(key=lambda item: (item.mismatches, item.start))
    deduped: List[MotifMatch] = []
    seen_windows = set()
    for match in matches:
        window_key = (match.start, match.end)
        if window_key in seen_windows:
            continue
        seen_windows.add(window_key)
        deduped.append(match)
    return deduped


def make_fragment_candidates(record: SequenceRecord) -> List[FragmentCandidate]:
    seq = record.cleaned_sequence
    if not seq:
        return []
    candidates: List[FragmentCandidate] = [
        FragmentCandidate(
            source="full_query",
            sequence=seq,
            insertion_anchor="left",
            note=(
                "Full cleaned query used directly; when no Tn5 motif is present, "
                "the query's first base is treated as the insertion-anchored end."
            ),
            source_priority=0,
        )
    ]
    for match in motif_matches(seq)[:3]:
        left = seq[: match.start]
        right = seq[match.end :]
        motif_note = f"{match.label}@{match.start + 1}-{match.end},mm={match.mismatches}"
        if len(left) >= MIN_FRAGMENT_LENGTH:
            candidates.append(
                FragmentCandidate(
                    source="left_of_tn5",
                    sequence=left,
                    insertion_anchor="right",
                    note=f"Left flank before {motif_note}. Insertion inferred from fragment end.",
                    motif_match=match,
                    source_priority=1,
                )
            )
        if len(right) >= MIN_FRAGMENT_LENGTH:
            candidates.append(
                FragmentCandidate(
                    source="right_of_tn5",
                    sequence=right,
                    insertion_anchor="left",
                    note=f"Right flank after {motif_note}. Insertion inferred from fragment start.",
                    motif_match=match,
                    source_priority=1,
                )
            )
    unique: Dict[Tuple[str, str, str], FragmentCandidate] = {}
    for candidate in candidates:
        key = (candidate.source, candidate.sequence, candidate.insertion_anchor)
        if key not in unique or candidate.note < unique[key].note:
            unique[key] = candidate
    return list(unique.values())


def seed_offsets(seq_len: int, k: int) -> List[int]:
    if seq_len < k:
        return []
    if seq_len == k:
        return [0]
    offsets = {
        0,
        max(0, seq_len // 4 - k // 2),
        max(0, seq_len // 2 - k // 2),
        max(0, (3 * seq_len) // 4 - k // 2),
        max(0, seq_len - k),
    }
    return sorted(offsets)


def identity_threshold(seq_len: int) -> float:
    return HIGH_IDENTITY_SHORT if seq_len < 60 else HIGH_IDENTITY_LONG


def rescue_identity_threshold(seq_len: int) -> float:
    return RESCUE_IDENTITY_SHORT if seq_len < 45 else RESCUE_IDENTITY_LONG


def score_hit(aligned_query: str, ref_seq: str) -> Tuple[int, int, int]:
    comparable_bases = 0
    matched_bases = 0
    for q_base, r_base in zip(aligned_query, ref_seq):
        if q_base not in "ACGT":
            continue
        comparable_bases += 1
        if q_base == r_base:
            matched_bases += 1
    return matched_bases, comparable_bases, len(aligned_query)


def junction_boundary_0based(
    strand: str,
    insertion_anchor: str,
    alignment_start: int,
    alignment_end: int,
) -> int:
    """Return the reference boundary between bases using 0-based coordinates.

    ``alignment_start``/``alignment_end`` are 1-based inclusive coordinates.
    The boundary is an index in [0, genome_length], i.e. the boundary before
    reference base boundary+1.  ``Insertion site`` remains the adjacent
    1-based reference base for backwards-compatible result tables.
    """
    start0 = alignment_start - 1
    end0 = alignment_end
    if strand == "+":
        return start0 if insertion_anchor == "left" else end0
    return end0 if insertion_anchor == "left" else start0


def _hit_from_ungapped_window(
    record: SequenceRecord,
    candidate: FragmentCandidate,
    reference_seqid: str,
    strand: str,
    aligned_query: str,
    start: int,
    genome: str,
) -> Optional[Hit]:
    ref_seq = genome[start : start + len(aligned_query)]
    matched_bases, comparable_bases, alignment_length = score_hit(aligned_query, ref_seq)
    if comparable_bases == 0:
        return None
    identity = matched_bases / comparable_bases
    query_coverage = comparable_bases / len(aligned_query)
    if identity < LOW_CONFIDENCE_IDENTITY or query_coverage < 0.50:
        return None
    alignment_start = start + 1
    alignment_end = start + len(aligned_query)
    if strand == "+":
        query_pos1_reference_coord = alignment_start
        insertion_site = alignment_start if candidate.insertion_anchor == "left" else alignment_end
    else:
        query_pos1_reference_coord = alignment_end
        insertion_site = alignment_end if candidate.insertion_anchor == "left" else alignment_start
    boundary = junction_boundary_0based(
        strand, candidate.insertion_anchor, alignment_start, alignment_end
    )
    note_parts = [candidate.note]
    if candidate.motif_match is None:
        note_parts.append("No Tn5 motif trimmed.")
    return Hit(
        query_id=record.query_id,
        row_number=record.row_number,
        candidate_source=candidate.source,
        reference_seqid=reference_seqid,
        strand=strand,
        alignment_start=alignment_start,
        alignment_end=alignment_end,
        query_pos1_reference_coord=query_pos1_reference_coord,
        insertion_site=insertion_site,
        identity=identity,
        query_coverage=query_coverage,
        alignment_length=alignment_length,
        matched_bases=matched_bases,
        comparable_bases=comparable_bases,
        source_priority=candidate.source_priority,
        note=" ".join(note_parts),
        alignment_score=float(2 * matched_bases - comparable_bases),
        edit_distance=comparable_bases - matched_bases,
        junction_boundary_0based=boundary,
        alignment_method="internal_exhaustive_ungapped",
    )


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def available_aligner(requested: str) -> str:
    if requested != "auto":
        if requested == "internal" or shutil.which(requested if requested != "blastn" else "blastn"):
            return requested
        raise RuntimeError(
            f"Requested aligner '{requested}' was not found on PATH. "
            "Install it or use --aligner internal."
        )
    for executable in ("minimap2", "blastn", "bwa"):
        if shutil.which(executable):
            return "blastn" if executable == "blastn" else executable
    return "internal"


def _sam_tag(fields: Sequence[str], tag_name: str) -> Optional[str]:
    prefix = tag_name + ":"
    for field in fields[11:]:
        if field.startswith(prefix):
            return field.split(":", 2)[-1]
    return None


def _external_hit_from_sam(
    record: SequenceRecord,
    candidate: FragmentCandidate,
    reference_seqid: str,
    aligner: str,
    fields: Sequence[str],
) -> Optional[Hit]:
    if len(fields) < 11 or fields[2] == "*":
        return None
    qname = fields[0]
    flag = int(fields[1])
    cigar = fields[5]
    if cigar == "*":
        return None
    ops = [(int(length), op) for length, op in _CIGAR_RE.findall(cigar)]
    q_aligned = sum(length for length, op in ops if op in "MI=X")
    q_start_clip = 0
    for length, op in ops:
        if op == "S":
            q_start_clip += length
        else:
            break
    qlen = len(candidate.sequence)
    if q_aligned == 0 or qlen == 0:
        return None
    q_start = q_start_clip + 1
    q_end = q_start_clip + q_aligned
    # An insertion coordinate is only supported when the aligned segment
    # reaches the insertion-anchored end of the query.
    if candidate.insertion_anchor == "left" and q_start != 1:
        return None
    if candidate.insertion_anchor == "right" and q_end != qlen:
        return None
    target_span = sum(length for length, op in ops if op in "MDN=X")
    alignment_start = int(fields[3])
    alignment_end = alignment_start + target_span - 1
    nm = int(_sam_tag(fields, "NM") or "0")
    insertions = sum(length for length, op in ops if op == "I")
    deletions = sum(length for length, op in ops if op in "DN")
    mismatches = max(0, nm - insertions - deletions)
    comparable_bases = q_aligned
    matched_bases = max(0, q_aligned - insertions - mismatches)
    identity = matched_bases / max(1, comparable_bases)
    query_coverage = q_aligned / qlen
    if identity < LOW_CONFIDENCE_IDENTITY or query_coverage < 0.50:
        return None
    plus_query = qname == "plus"
    reverse_ref = bool(flag & 16)
    strand = "+" if (plus_query == (not reverse_ref)) else "-"
    if strand == "+":
        query_pos1_reference_coord = alignment_start
        insertion_site = alignment_start if candidate.insertion_anchor == "left" else alignment_end
    else:
        query_pos1_reference_coord = alignment_end
        insertion_site = alignment_end if candidate.insertion_anchor == "left" else alignment_start
    boundary = junction_boundary_0based(
        strand, candidate.insertion_anchor, alignment_start, alignment_end
    )
    score_tag = _sam_tag(fields, "AS")
    alignment_score = float(score_tag) if score_tag is not None else float(2 * matched_bases - comparable_bases)
    note_parts = [candidate.note, f"Standard local alignment via {aligner} SAM."]
    if candidate.motif_match is None:
        note_parts.append("No Tn5 motif trimmed.")
    return Hit(
        query_id=record.query_id,
        row_number=record.row_number,
        candidate_source=candidate.source,
        reference_seqid=reference_seqid,
        strand=strand,
        alignment_start=alignment_start,
        alignment_end=alignment_end,
        query_pos1_reference_coord=query_pos1_reference_coord,
        insertion_site=insertion_site,
        identity=identity,
        query_coverage=query_coverage,
        alignment_length=q_aligned,
        matched_bases=matched_bases,
        comparable_bases=comparable_bases,
        source_priority=candidate.source_priority,
        note=" ".join(note_parts),
        alignment_score=alignment_score,
        edit_distance=nm,
        junction_boundary_0based=boundary,
        alignment_method=aligner,
    )


def _external_hit_from_blast(
    record: SequenceRecord,
    candidate: FragmentCandidate,
    reference_seqid: str,
    fields: Sequence[str],
) -> Optional[Hit]:
    # qseqid sseqid sstart send length nident qlen qstart qend bitscore evalue
    if len(fields) < 11:
        return None
    qname, _, sstart, send, aln_len, nident, qlen, qstart, qend, bitscore, _ = fields[:11]
    qlen_i = int(qlen)
    qstart_i, qend_i = int(qstart), int(qend)
    if candidate.insertion_anchor == "left" and qstart_i != 1:
        return None
    if candidate.insertion_anchor == "right" and qend_i != qlen_i:
        return None
    sstart_i, send_i = int(sstart), int(send)
    plus_query = qname == "plus"
    reference_reverse = sstart_i > send_i
    strand = "+" if (plus_query == (not reference_reverse)) else "-"
    alignment_start, alignment_end = min(sstart_i, send_i), max(sstart_i, send_i)
    identity = int(nident) / max(1, int(aln_len))
    query_coverage = (qend_i - qstart_i + 1) / max(1, qlen_i)
    if identity < LOW_CONFIDENCE_IDENTITY or query_coverage < 0.50:
        return None
    if strand == "+":
        query_pos1_reference_coord = alignment_start
        insertion_site = alignment_start if candidate.insertion_anchor == "left" else alignment_end
    else:
        query_pos1_reference_coord = alignment_end
        insertion_site = alignment_end if candidate.insertion_anchor == "left" else alignment_start
    boundary = junction_boundary_0based(
        strand, candidate.insertion_anchor, alignment_start, alignment_end
    )
    comparable_bases = int(aln_len)
    matched_bases = int(nident)
    note_parts = [candidate.note, "Standard local alignment via BLASTn-short."]
    if candidate.motif_match is None:
        note_parts.append("No Tn5 motif trimmed.")
    return Hit(
        query_id=record.query_id,
        row_number=record.row_number,
        candidate_source=candidate.source,
        reference_seqid=reference_seqid,
        strand=strand,
        alignment_start=alignment_start,
        alignment_end=alignment_end,
        query_pos1_reference_coord=query_pos1_reference_coord,
        insertion_site=insertion_site,
        identity=identity,
        query_coverage=query_coverage,
        alignment_length=int(aln_len),
        matched_bases=matched_bases,
        comparable_bases=comparable_bases,
        source_priority=candidate.source_priority,
        note=" ".join(note_parts),
        alignment_score=float(bitscore),
        edit_distance=comparable_bases - matched_bases,
        junction_boundary_0based=boundary,
        alignment_method="blastn-short",
    )


def ensure_bwa_index(genome_path: Path) -> None:
    """Build the classic BWA index once when it is not already present."""
    suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
    index_paths = [Path(str(genome_path) + suffix) for suffix in suffixes]
    if all(path.is_file() and path.stat().st_size > 0 for path in index_paths):
        return
    try:
        subprocess.run(
            ["bwa", "index", str(genome_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            f"BWA index construction failed for {genome_path}: {exc}"
        ) from exc


def map_fragment_external(
    record: SequenceRecord,
    candidate: FragmentCandidate,
    reference_seqid: str,
    genome_path: Path,
    aligner: str,
) -> List[Hit]:
    """Run one standard local aligner against both query orientations."""
    with tempfile.TemporaryDirectory(prefix="tn5-map-") as temp_dir:
        query_path = Path(temp_dir) / "query.fa"
        query_path.write_text(
            f">plus\n{candidate.sequence}\n>minus\n{reverse_complement(candidate.sequence)}\n",
            encoding="ascii",
        )
        if aligner == "blastn":
            command = [
                "blastn",
                "-task",
                "blastn-short",
                "-query",
                str(query_path),
                "-subject",
                str(genome_path),
                "-dust",
                "no",
                "-outfmt",
                "6 qseqid sseqid sstart send length nident qlen qstart qend bitscore evalue",
            ]
        elif aligner == "minimap2":
            command = [
                "minimap2",
                "-a",
                "--MD",
                "-N",
                "50",
                "-x",
                "sr",
                str(genome_path),
                str(query_path),
            ]
        elif aligner == "bwa":
            ensure_bwa_index(genome_path)
            command = ["bwa", "mem", "-a", str(genome_path), str(query_path)]
        else:
            raise ValueError(f"Unsupported external aligner: {aligner}")
        try:
            completed = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            raise RuntimeError(
                f"{aligner} failed while mapping {record.query_id}: {exc}"
            ) from exc
    hits_by_key: Dict[Tuple[str, int], Hit] = {}
    for line in completed.stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t") if aligner != "blastn" else line.split("\t")
        hit = (
            _external_hit_from_blast(record, candidate, reference_seqid, fields)
            if aligner == "blastn"
            else _external_hit_from_sam(
                record,
                candidate,
                reference_seqid,
                aligner,
                fields,
            )
        )
        if hit is None:
            continue
        key = (hit.strand, hit.alignment_start)
        previous = hits_by_key.get(key)
        if previous is None or (hit.alignment_score, hit.identity, hit.query_coverage) > (
            previous.alignment_score,
            previous.identity,
            previous.query_coverage,
        ):
            hits_by_key[key] = hit
    hits = sorted(
        hits_by_key.values(),
        key=lambda item: (item.alignment_score, item.identity, item.query_coverage),
        reverse=True,
    )
    top_score = hits[0].alignment_score if hits else 0.0
    second_score = hits[1].alignment_score if len(hits) > 1 else None
    for rank, hit in enumerate(hits, start=1):
        hit.hit_rank = rank
        hit.score_gap = (
            top_score - second_score
            if rank == 1 and second_score is not None
            else top_score - hit.alignment_score
        )
    return hits


def map_fragment(
    record: SequenceRecord,
    candidate: FragmentCandidate,
    reference_seqid: str,
    genome: str,
    indices: Dict[int, Dict[str, List[int]]],
    seed_sizes: Optional[Sequence[int]] = None,
    aligner: str = "internal",
    genome_path: Optional[Path] = None,
    alignment_cache: Optional[Dict[Tuple[str, str, str], List[Hit]]] = None,
) -> List[Hit]:
    seq = candidate.sequence
    if len(seq) < MIN_FRAGMENT_LENGTH:
        return []
    if aligner != "internal" and genome_path is not None:
        cache_key = (aligner, candidate.sequence, candidate.insertion_anchor)
        if alignment_cache is not None and cache_key in alignment_cache:
            cached_hits = alignment_cache[cache_key]
            return [
                replace(
                    hit,
                    query_id=record.query_id,
                    row_number=record.row_number,
                    candidate_source=candidate.source,
                    source_priority=candidate.source_priority,
                    note=hit.note.replace(hit.note.split(";", 1)[0], candidate.note, 1),
                )
                for hit in cached_hits
            ]
        external_hits = map_fragment_external(
            record, candidate, reference_seqid, genome_path, aligner
        )
        if alignment_cache is not None:
            alignment_cache[cache_key] = [
                replace(hit, query_id="", row_number=0) for hit in external_hits
            ]
        return external_hits
    active_seed_sizes = tuple(seed_sizes or sorted(indices.keys(), reverse=True))
    hits_by_key: Dict[Tuple[str, int], Hit] = {}
    for strand, aligned_query in (("+", seq), ("-", reverse_complement(seq))):
        all_candidate_positions = set()
        for k in active_seed_sizes:
            if len(aligned_query) < k:
                continue
            index = indices[k]
            for offset in seed_offsets(len(aligned_query), k):
                seed = aligned_query[offset : offset + k]
                if "N" in seed:
                    continue
                positions = index.get(seed)
                # Very high-copy seeds (e.g. rRNA/simple repeats) are not
                # informative and can generate millions of windows. They are
                # skipped, but every informative seed size is still searched;
                # importantly, we no longer stop after the first seed size.
                if not positions or len(positions) > 2500:
                    continue
                for pos in positions:
                    start = pos - offset
                    if 0 <= start <= len(genome) - len(aligned_query):
                        all_candidate_positions.add(start)
        for start in all_candidate_positions:
            hit = _hit_from_ungapped_window(
                record, candidate, reference_seqid, strand, aligned_query, start, genome
            )
            if hit is None:
                continue
            key = (strand, start)
            previous = hits_by_key.get(key)
            if previous is None or (
                hit.alignment_score,
                hit.identity,
                hit.query_coverage,
                candidate.source_priority,
            ) > (
                previous.alignment_score,
                previous.identity,
                previous.query_coverage,
                previous.source_priority,
            ):
                hits_by_key[key] = hit
    hits = sorted(
        hits_by_key.values(),
        key=lambda item: (
            item.alignment_score,
            item.matched_bases,
            item.identity,
            item.query_coverage,
            item.alignment_length,
        ),
        reverse=True,
    )
    top_score = hits[0].alignment_score if hits else 0.0
    second_score = hits[1].alignment_score if len(hits) > 1 else None
    for rank, hit in enumerate(hits, start=1):
        hit.hit_rank = rank
        hit.score_gap = (
            top_score - second_score
            if rank == 1 and second_score is not None
            else top_score - hit.alignment_score
        )
    return hits


def rescue_subfragments(candidate: FragmentCandidate) -> List[FragmentCandidate]:
    seq = candidate.sequence
    fragments: List[FragmentCandidate] = []
    for window_size in RESCUE_WINDOW_SIZES:
        span = min(window_size, len(seq))
        if span < RESCUE_MIN_SPAN:
            continue
        if candidate.insertion_anchor == "left":
            sub_seq = seq[:span]
            side_label = "prefix"
        else:
            sub_seq = seq[-span:]
            side_label = "suffix"
        fragments.append(
            FragmentCandidate(
                source=f"rescue_{side_label}_{span}_from_{candidate.source}",
                sequence=sub_seq,
                insertion_anchor=candidate.insertion_anchor,
                note=(
                    f"Rescue mode used the insertion-anchored {side_label} subfragment "
                    f"({span} bp) derived from {candidate.source}."
                ),
                motif_match=candidate.motif_match,
                source_priority=candidate.source_priority + 10,
            )
        )
        if span == len(seq):
            break
    unique: Dict[str, FragmentCandidate] = {}
    for fragment in fragments:
        unique[fragment.source] = fragment
    return list(unique.values())


def evaluate_hits(candidate: FragmentCandidate, hits: List[Hit]) -> CandidateEvaluation:
    if not hits:
        return CandidateEvaluation(
            candidate=candidate,
            hits=[],
            status="no_hit",
            top_hit=None,
            note="No candidate alignment passed the minimum identity cutoff.",
        )
    top = hits[0]
    top_cutoff = identity_threshold(len(candidate.sequence))
    equal_best = [
        hit
        for hit in hits
        if hit.matched_bases == top.matched_bases and hit.comparable_bases == top.comparable_bases
    ]
    if top.identity >= top_cutoff and top.query_coverage >= MIN_QUERY_COVERAGE:
        if len(equal_best) > 1:
            for hit in hits:
                hit.hit_status = "multiple_hits"
                hit.hit_count = len(equal_best)
            return CandidateEvaluation(
                candidate=candidate,
                hits=hits,
                status="multiple_hits",
                top_hit=top,
                note="Multiple equally scoring placements remained after conservative filtering.",
            )
        if len(hits) > 1 and (
            (top.matched_bases - hits[1].matched_bases) < TOP_HIT_MARGIN
            or (top.alignment_score - hits[1].alignment_score) < TOP_HIT_SCORE_GAP
        ):
            for hit in hits:
                hit.hit_status = "multiple_hits"
                hit.hit_count = len(hits)
            return CandidateEvaluation(
                candidate=candidate,
                hits=hits,
                status="multiple_hits",
                top_hit=top,
                note="The top two placements differed by fewer than two matched bases.",
            )
        for hit in hits:
            hit.hit_status = "unique_hit" if hit.hit_rank == 1 else "secondary_hit"
            hit.hit_count = len(hits)
        return CandidateEvaluation(
            candidate=candidate,
            hits=hits,
            status="unique_hit",
            top_hit=top,
            note="A single high-confidence placement was retained.",
        )
    for hit in hits:
        hit.hit_status = "low_confidence"
        hit.hit_count = len(hits)
    return CandidateEvaluation(
        candidate=candidate,
        hits=hits,
        status="low_confidence",
        top_hit=top,
        note="Best placement did not meet the conservative identity and coverage cutoffs.",
    )


def evaluate_rescue_hits(candidate: FragmentCandidate, hits: List[Hit]) -> CandidateEvaluation:
    if not hits:
        return CandidateEvaluation(
            candidate=candidate,
            hits=[],
            status="no_hit",
            top_hit=None,
            note="No rescue placement passed the minimum cutoff.",
        )
    top = hits[0]
    top_cutoff = rescue_identity_threshold(len(candidate.sequence))
    equal_best = [
        hit
        for hit in hits
        if hit.matched_bases == top.matched_bases and hit.comparable_bases == top.comparable_bases
    ]
    if top.identity >= top_cutoff and top.query_coverage >= 1.0:
        if len(equal_best) > 1:
            for hit in hits:
                hit.hit_status = "multiple_hits"
                hit.hit_count = len(equal_best)
            return CandidateEvaluation(
                candidate=candidate,
                hits=hits,
                status="multiple_hits",
                top_hit=top,
                note="Rescue mode found multiple equally scoring insertion-anchored placements.",
            )
        if len(hits) > 1 and (
            (top.matched_bases - hits[1].matched_bases) < TOP_HIT_MARGIN
            or (top.alignment_score - hits[1].alignment_score) < TOP_HIT_SCORE_GAP
        ):
            for hit in hits:
                hit.hit_status = "multiple_hits"
                hit.hit_count = len(hits)
            return CandidateEvaluation(
                candidate=candidate,
                hits=hits,
                status="multiple_hits",
                top_hit=top,
                note="Rescue mode top placements remained too close to separate confidently.",
            )
        for hit in hits:
            hit.hit_status = "rescued_unique_hit" if hit.hit_rank == 1 else "secondary_hit"
            hit.hit_count = len(hits)
        return CandidateEvaluation(
            candidate=candidate,
            hits=hits,
            status="rescued_unique_hit",
            top_hit=top,
            note="Insertion site was rescued by a unique insertion-anchored subfragment placement.",
        )
    for hit in hits:
        hit.hit_status = "low_confidence"
        hit.hit_count = len(hits)
    return CandidateEvaluation(
        candidate=candidate,
        hits=hits,
        status="low_confidence",
        top_hit=top,
        note="Rescue mode still failed the identity or uniqueness cutoff.",
    )


def choose_consistent_rescue_evaluation(
    evaluations: List[CandidateEvaluation],
) -> Optional[CandidateEvaluation]:
    """Require all qualifying rescue windows to agree on one junction."""
    qualified = [item for item in evaluations if item.status == "rescued_unique_hit" and item.top_hit]
    if not qualified:
        return None
    groups: Dict[Tuple[str, str, int], List[CandidateEvaluation]] = defaultdict(list)
    for item in qualified:
        hit = item.top_hit
        assert hit is not None
        groups[(hit.reference_seqid, hit.strand, hit.junction_boundary_0based)].append(item)
    if len(groups) > 1:
        best = max(
            qualified,
            key=lambda item: (
                item.top_hit.alignment_score if item.top_hit else -1.0,
                item.top_hit.identity if item.top_hit else -1.0,
                len(item.candidate.sequence),
            ),
        )
        conflict_hits = [item.top_hit for item in qualified if item.top_hit is not None]
        for hit in conflict_hits:
            hit.hit_status = "multiple_hits"
            hit.hit_count = len(conflict_hits)
        best.hits = conflict_hits
        best.top_hit = max(
            conflict_hits,
            key=lambda hit: (hit.alignment_score, hit.identity, hit.query_coverage),
        )
        best.status = "multiple_hits"
        best.note = (
            "Multiple qualifying rescue windows mapped to conflicting junction "
            "coordinates or orientations; rescue was marked ambiguous."
        )
        return best
    return max(
        qualified,
        key=lambda item: (
            item.top_hit.alignment_score if item.top_hit else -1.0,
            item.top_hit.identity if item.top_hit else -1.0,
            len(item.candidate.sequence),
        ),
    )


def choose_best_evaluation(evaluations: List[CandidateEvaluation]) -> CandidateEvaluation:
    status_rank = {
        "rescued_unique_hit": 4,
        "unique_hit": 3,
        "multiple_hits": 2,
        "low_confidence": 1,
        "no_hit": 0,
    }
    return max(
        evaluations,
        key=lambda item: (
            status_rank[item.status],
            item.top_hit.matched_bases if item.top_hit else -1,
            item.top_hit.identity if item.top_hit else -1.0,
            len(item.candidate.sequence),
            item.candidate.source_priority,
        ),
    )


def overlapping_features(
    features: Sequence[GeneFeature] | Sequence[CDSFeature], coord: int
) -> List[GeneFeature] | List[CDSFeature]:
    return [feature for feature in features if feature.start <= coord <= feature.end]


def nearest_genes(genes: Sequence[GeneFeature], coord: int) -> Tuple[Optional[GeneFeature], Optional[GeneFeature]]:
    upstream = None
    downstream = None
    for gene in genes:
        if gene.end < coord and (upstream is None or gene.end > upstream.end):
            upstream = gene
        if gene.start > coord and (downstream is None or gene.start < downstream.start):
            downstream = gene
    return upstream, downstream


def outside_gene_relation(gene: GeneFeature, coord: int) -> Tuple[str, int, str]:
    """Describe an external position relative to a gene's transcription direction."""
    if coord < gene.start:
        genomic_side = "left"
        distance = gene.start - coord
    elif coord > gene.end:
        genomic_side = "right"
        distance = coord - gene.end
    else:
        return "within_gene", 0, "overlap"

    if gene.strand == "+":
        relation = "5prime_upstream" if genomic_side == "left" else "3prime_downstream"
    elif gene.strand == "-":
        relation = "3prime_downstream" if genomic_side == "left" else "5prime_upstream"
    else:
        relation = f"unknown_orientation_{genomic_side}"
    return relation, distance, genomic_side


def intergenic_distance_score(relation: str, distance: int) -> Tuple[int, str]:
    """Score a reported candidate; 500 bp remains only the outer reporting limit."""
    if relation != "5prime_upstream":
        return 0, "3prime_context_only"
    if distance <= 100:
        return 40, "5prime_distance_0_100"
    if distance <= 300:
        return 25, "5prime_distance_101_300"
    return 10, "5prime_distance_301_500"


def score_priority(score: int) -> str:
    if score >= 40:
        return "high"
    if score >= 25:
        return "moderate"
    if score >= 10:
        return "low"
    return "context_only"


def profile_has_outward_promoters(profile: str) -> bool:
    return profile in {
        "bidirectional_outward_promoters",
        "bidirectional_promoters_and_terminator",
    }


def profile_has_terminator(profile: str) -> bool:
    return profile in {
        "transcriptional_terminator",
        "bidirectional_promoters_and_terminator",
    }


def cds_for_gene(gene: GeneFeature, cdss: Sequence[CDSFeature]) -> Optional[CDSFeature]:
    identifiers = {gene.locus_tag, gene.gene_id, gene.gene_name}
    identifiers.discard("")
    for cds in cdss:
        cds_identifiers = {cds.locus_tag, cds.gene_id, cds.gene_name}
        cds_identifiers.discard("")
        if identifiers & cds_identifiers:
            return cds
    return None


def candidate_fields(
    candidates: Sequence[NearbyGeneCandidate],
    cdss: Sequence[CDSFeature],
) -> Dict[str, str]:
    separator = " | "
    names: List[str] = []
    codes: List[str] = []
    relations: List[str] = []
    distances: List[str] = []
    strands: List[str] = []
    priorities: List[str] = []
    functions: List[str] = []
    ontology_terms: List[str] = []
    go_tags: List[str] = []
    scores: List[str] = []
    evidence: List[str] = []
    operon_ids: List[str] = []
    operon_relations: List[str] = []
    for candidate in candidates:
        gene = candidate.gene
        cds = cds_for_gene(gene, cdss)
        names.append(gene.gene_name)
        codes.append(gene.gene_code)
        relations.append(candidate.relation)
        distances.append(str(candidate.distance))
        strands.append(gene.strand)
        priorities.append(candidate.priority)
        functions.append(cds.product if cds else "")
        ontology_terms.append(join_unique(cds.ontology_terms) if cds else "")
        go_tags.append(join_unique(cds.go_tags) if cds else "")
        scores.append(str(candidate.score))
        evidence.append("+".join(candidate.evidence))
        operon_ids.append(candidate.operon_id or "NA")
        operon_relations.append(candidate.operon_relation)
    return {
        "candidate_gene_name": separator.join(names),
        "candidate_gene_code": separator.join(codes),
        "candidate_relation": separator.join(relations),
        "candidate_distance_bp": separator.join(distances),
        "candidate_strand": separator.join(strands),
        "candidate_priority": separator.join(priorities),
        "candidate_function": separator.join(functions),
        # GO annotations contain " | " between a term name and GO identifier.
        # Use a distinct candidate separator so downstream selection can preserve
        # the one-to-one candidate-to-annotation relationship.
        "candidate_ontology_term": " || ".join(ontology_terms),
        "candidate_go_tags": " || ".join(go_tags),
        "candidate_score": separator.join(scores),
        "candidate_evidence": separator.join(evidence),
        "candidate_operon_id": separator.join(operon_ids),
        "candidate_operon_relation": separator.join(operon_relations),
    }


def build_intergenic_candidates(
    flanking_genes: Sequence[GeneFeature],
    coord: int,
    window: int,
    operon_map: Dict[str, str],
    tn5_structure_profile: str,
) -> List[NearbyGeneCandidate]:
    if window <= 0:
        return []

    in_window: List[Tuple[GeneFeature, str, int]] = []
    for gene in flanking_genes:
        relation, distance, _ = outside_gene_relation(gene, coord)
        if distance <= window:
            in_window.append((gene, relation, distance))

    operon_ids = {id(gene): operon_id_for_gene(gene, operon_map) for gene in flanking_genes}
    if not operon_map:
        shared_operon = ""
        flank_operon_status = "not_assessed"
    elif len(flanking_genes) < 2 or any(not operon_ids[id(gene)] for gene in flanking_genes):
        shared_operon = ""
        flank_operon_status = "insufficient_annotation"
    elif len({gene.strand for gene in flanking_genes}) != 1:
        shared_operon = ""
        flank_operon_status = "different_transcription_strands"
    elif len({operon_ids[id(gene)] for gene in flanking_genes}) == 1:
        shared_operon = operon_ids[id(flanking_genes[0])]
        flank_operon_status = "same_operon"
    else:
        shared_operon = ""
        flank_operon_status = "different_operons"

    candidates: List[NearbyGeneCandidate] = []
    for gene, relation, distance in in_window:
        score, distance_evidence = intergenic_distance_score(relation, distance)
        evidence = [distance_evidence]
        operon_id = operon_ids.get(id(gene), "")
        operon_relation = flank_operon_status

        if shared_operon:
            if relation == "5prime_upstream":
                operon_relation = "same_operon_downstream"
                score += 20
                evidence.append("same_operon_polar_candidate")
            else:
                operon_relation = "same_operon_upstream_context"
                evidence.append("same_operon_upstream_context")

        if profile_has_outward_promoters(tn5_structure_profile) and relation == "5prime_upstream":
            score += 10
            evidence.append("bidirectional_outward_promoter_toward_gene")
        if (
            profile_has_terminator(tn5_structure_profile)
            and operon_relation == "same_operon_downstream"
        ):
            score += 10
            evidence.append("terminator_polar_effect")
        if tn5_structure_profile == "unknown":
            evidence.append("tn5_structure_unknown")

        candidates.append(
            NearbyGeneCandidate(
                gene=gene,
                relation=relation,
                distance=distance,
                score=score,
                priority=score_priority(score),
                evidence=evidence,
                operon_id=operon_id,
                operon_relation=operon_relation,
            )
        )

    candidates.sort(
        key=lambda item: (
            -item.score,
            item.distance,
            item.gene.gene_code or item.gene.gene_name,
        )
    )
    return candidates


def join_unique(values: Iterable[str]) -> str:
    seen = set()
    ordered = []
    for value in values:
        if not value or value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return "; ".join(ordered)


def annotate_insertion(
    seqid: str,
    coord: int,
    genes: Dict[str, List[GeneFeature]],
    cdss: Dict[str, List[CDSFeature]],
    gene_lookup: Dict[Tuple[str, str], GeneFeature],
    intergenic_candidate_window: int = 500,
    operon_map: Optional[Dict[str, str]] = None,
    tn5_structure_profile: str = "unknown",
) -> AnnotationResult:
    operon_map = operon_map or {}
    seq_genes = genes.get(seqid, [])
    seq_cdss = cdss.get(seqid, [])
    overlapping_cdss = list(overlapping_features(seq_cdss, coord))
    overlapping_genes = list(overlapping_features(seq_genes, coord))
    if overlapping_cdss:
        cds = overlapping_cdss[0]
        gene = (
            gene_lookup.get((seqid, cds.locus_tag))
            or gene_lookup.get((seqid, cds.gene_id))
            or gene_lookup.get((seqid, cds.gene_name))
        )
        gene_name = gene.gene_name if gene else cds.gene_name
        gene_code = gene.gene_code if gene else (cds.locus_tag or cds.gene_id)
        note = ""
        if len(overlapping_cdss) > 1:
            note = "Multiple CDS features overlap this insertion; the first CDS in coordinate order was reported."
        return AnnotationResult(
            feature_relation="CDS",
            gene_name=gene_name,
            gene_code=gene_code,
            putative_function=cds.product,
            ontology_term=join_unique(cds.ontology_terms),
            go_tags=join_unique(cds.go_tags),
            note=note,
            candidate_gene_name=gene_name,
            candidate_gene_code=gene_code,
            candidate_relation="within_CDS",
            candidate_distance_bp="0",
            candidate_strand=cds.strand,
            candidate_priority="primary",
            candidate_function=cds.product,
            candidate_ontology_term=join_unique(cds.ontology_terms),
            candidate_go_tags=join_unique(cds.go_tags),
            candidate_score="100",
            candidate_evidence="direct_CDS_overlap",
            candidate_operon_id=operon_id_for_gene(gene, operon_map) if gene else "",
            candidate_operon_relation="direct_gene",
            tn5_structure_profile=tn5_structure_profile,
        )
    if overlapping_genes:
        gene = overlapping_genes[0]
        gene_cds = cds_for_gene(gene, seq_cdss)
        note = "Insertion falls within a gene interval but outside annotated CDS features."
        if len(overlapping_genes) > 1:
            note += " Multiple gene intervals overlap this position; the first gene in coordinate order was reported."
        return AnnotationResult(
            feature_relation="gene_non_CDS",
            gene_name=gene.gene_name,
            gene_code=gene.gene_code,
            putative_function=gene_cds.product if gene_cds else "",
            ontology_term=join_unique(gene_cds.ontology_terms) if gene_cds else "",
            go_tags=join_unique(gene_cds.go_tags) if gene_cds else "",
            note=note,
            candidate_gene_name=gene.gene_name,
            candidate_gene_code=gene.gene_code,
            candidate_relation="within_gene_non_CDS",
            candidate_distance_bp="0",
            candidate_strand=gene.strand,
            candidate_priority="primary",
            candidate_function=gene_cds.product if gene_cds else "",
            candidate_ontology_term=(
                join_unique(gene_cds.ontology_terms) if gene_cds else ""
            ),
            candidate_go_tags=join_unique(gene_cds.go_tags) if gene_cds else "",
            candidate_score="90",
            candidate_evidence="direct_gene_interval_overlap",
            candidate_operon_id=operon_id_for_gene(gene, operon_map),
            candidate_operon_relation="direct_gene",
            tn5_structure_profile=tn5_structure_profile,
        )
    upstream, downstream = nearest_genes(seq_genes, coord)
    note_parts = ["Intergenic insertion."]
    flanking_genes = [gene for gene in (upstream, downstream) if gene is not None]
    if upstream is not None:
        relation, distance, _ = outside_gene_relation(upstream, coord)
        note_parts.append(
            f"Upstream gene {upstream.gene_name or upstream.gene_code} ends at {upstream.end} "
            f"({distance} bp away; {relation} relative to its {upstream.strand}-strand transcription)."
        )
    if downstream is not None:
        relation, distance, _ = outside_gene_relation(downstream, coord)
        note_parts.append(
            f"Downstream gene {downstream.gene_name or downstream.gene_code} starts at {downstream.start} "
            f"({distance} bp away; {relation} relative to its {downstream.strand}-strand transcription)."
        )
    nearby_candidates = build_intergenic_candidates(
        flanking_genes,
        coord,
        intergenic_candidate_window,
        operon_map,
        tn5_structure_profile,
    )
    fields = candidate_fields(nearby_candidates, seq_cdss)
    if nearby_candidates:
        note_parts.append(
            f"Genes within the {intergenic_candidate_window}-bp maximum reporting window "
            "were ranked by transcriptional side, distance tier, available operon "
            f"annotation, and Tn5 structure profile ({tn5_structure_profile}). "
            "The window is not an effect threshold, and candidate rank does not establish causality."
        )
    elif intergenic_candidate_window > 0:
        note_parts.append(
            f"No flanking gene falls within the {intergenic_candidate_window}-bp "
            "maximum candidate reporting window."
        )
    return AnnotationResult(
        feature_relation="intergenic",
        gene_name="",
        gene_code="",
        putative_function="",
        ontology_term="",
        go_tags="",
        note=" ".join(note_parts),
        tn5_structure_profile=tn5_structure_profile,
        **fields,
    )


def xml_col_name(index: int) -> str:
    result = []
    current = index
    while current:
        current, remainder = divmod(current - 1, 26)
        result.append(chr(65 + remainder))
    return "".join(reversed(result))


def worksheet_xml(rows: List[List[object]]) -> str:
    xml_rows = []
    for row_idx, row in enumerate(rows, start=1):
        cells = []
        for col_idx, value in enumerate(row, start=1):
            cell_ref = f"{xml_col_name(col_idx)}{row_idx}"
            if value is None:
                continue
            if isinstance(value, float):
                cell_xml = f'<c r="{cell_ref}"><v>{value:.6f}</v></c>'
            elif isinstance(value, int):
                cell_xml = f'<c r="{cell_ref}"><v>{value}</v></c>'
            else:
                text = escape(str(value))
                cell_xml = f'<c r="{cell_ref}" t="inlineStr"><is><t>{text}</t></is></c>'
            cells.append(cell_xml)
        xml_rows.append(f'<row r="{row_idx}">{"".join(cells)}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        "<sheetData>"
        + "".join(xml_rows)
        + "</sheetData></worksheet>"
    )


def write_xlsx(path: Path, sheets: List[Tuple[str, List[List[object]]]]) -> None:
    workbook_xml = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
        "<sheets>",
    ]
    workbook_rels = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    ]
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    ]
    for idx, (name, _) in enumerate(sheets, start=1):
        workbook_xml.append(f'<sheet name="{escape(name)}" sheetId="{idx}" r:id="rId{idx}"/>')
        workbook_rels.append(
            f'<Relationship Id="rId{idx}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{idx}.xml"/>'
        )
        content_types.append(
            f'<Override PartName="/xl/worksheets/sheet{idx}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
    workbook_xml.append("</sheets></workbook>")
    workbook_rels.append(
        '<Relationship Id="rId_styles" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
    )
    workbook_rels.append("</Relationships>")
    content_types.append("</Types>")
    root_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
        '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
        '<borders count="1"><border/></borders>'
        '<cellStyleXfs count="1"><xf/></cellStyleXfs>'
        '<cellXfs count="1"><xf xfId="0"/></cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        "</styleSheet>"
    )
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", "".join(content_types))
        zf.writestr("_rels/.rels", root_rels)
        zf.writestr("xl/workbook.xml", "".join(workbook_xml))
        zf.writestr("xl/_rels/workbook.xml.rels", "".join(workbook_rels))
        zf.writestr("xl/styles.xml", styles_xml)
        for idx, (_, rows) in enumerate(sheets, start=1):
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", worksheet_xml(rows))


def rows_from_dicts(headers: List[str], records: List[Dict[str, object]]) -> List[List[object]]:
    rows: List[List[object]] = [headers]
    for record in records:
        rows.append([record.get(header, "") for header in headers])
    return rows


def write_tsv(path: Path, headers: List[str], rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    if not 0 <= args.intergenic_candidate_window <= 500:
        raise ValueError(
            "--intergenic-candidate-window must be between 0 and 500 bp; "
            "500 bp is the maximum reporting window, not an effect threshold."
        )
    query_csv = Path(args.query_csv)
    gtf_path = Path(args.gtf)
    genome_path = Path(args.genome_fasta)
    operon_path = Path(args.operon_tsv) if args.operon_tsv else None
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / args.output_prefix

    reference_seqid, genome = read_single_fasta(genome_path)
    aligner = available_aligner(args.aligner)
    queries = load_queries(query_csv)
    genes, cdss = parse_gtf(gtf_path)
    gene_lookup = build_gene_lookup(genes)
    operon_map = load_operon_map(operon_path)
    indices = {k: build_kmer_index(genome, k) for k in SEED_SIZES}

    main_rows: List[Dict[str, object]] = []
    raw_match_rows: List[Dict[str, object]] = []
    failed_rows: List[Dict[str, object]] = []
    summary_counts = Counter()
    alignment_cache: Dict[Tuple[str, str, str], List[Hit]] = {}

    for record in queries:
        base_notes = list(record.notes)
        if not record.cleaned_sequence:
            summary_counts["total"] += 1
            summary_counts["no_sequence"] += 1
            row = {
                "Query_ID": record.query_id,
                "Row_number": record.row_number,
                "Gene name": "",
                "Gene code": "",
                "Insertion site": "",
                "Putative function or description": "",
                "Ontology_term": "",
                "Go_tags": "",
                "Feature_relation": "",
                "Reference_seqid": "",
                "Strand": "",
                "Hit_status": "no_sequence",
                "Note": join_unique(base_notes + ["No sequence remained after cleaning."]),
            }
            main_rows.append(row)
            failed_rows.append(
                {
                    "Query_ID": record.query_id,
                    "Row_number": record.row_number,
                    "Category": "no_sequence",
                    "Reason": row["Note"],
                }
            )
            continue

        summary_counts["total"] += 1
        summary_counts["valid_sequences"] += 1
        candidates = make_fragment_candidates(record)
        evaluations = [
            evaluate_hits(
                candidate,
                map_fragment(
                    record,
                    candidate,
                    reference_seqid,
                    genome,
                    indices,
                    SEED_SIZES,
                    aligner,
                    genome_path,
                    alignment_cache,
                ),
            )
            for candidate in candidates
        ]
        chosen = choose_best_evaluation(evaluations)
        rescued = False
        if chosen.status in {"low_confidence", "no_hit", "multiple_hits"}:
            rescue_candidates: List[FragmentCandidate] = []
            for candidate in candidates:
                rescue_candidates.extend(rescue_subfragments(candidate))
            rescue_indices = {k: indices.get(k) or build_kmer_index(genome, k) for k in RESCUE_SEED_SIZES}
            rescue_evaluations = [
                evaluate_rescue_hits(
                    rescue_candidate,
                    map_fragment(
                        record,
                        rescue_candidate,
                        reference_seqid,
                        genome,
                        rescue_indices,
                        RESCUE_SEED_SIZES,
                        aligner,
                        genome_path,
                        alignment_cache,
                    ),
                )
                for rescue_candidate in rescue_candidates
            ]
            if rescue_evaluations:
                rescued_choice = choose_consistent_rescue_evaluation(rescue_evaluations)
                if rescued_choice is not None and rescued_choice.status == "rescued_unique_hit":
                    chosen = rescued_choice
                    rescued = True
                elif rescued_choice is not None and rescued_choice.status == "multiple_hits":
                    chosen = rescued_choice
        chosen_note = join_unique(base_notes + [chosen.candidate.note, chosen.note])

        if chosen.status in {"unique_hit", "rescued_unique_hit"} and chosen.top_hit is not None:
            annotation = annotate_insertion(
                chosen.top_hit.reference_seqid,
                chosen.top_hit.insertion_site,
                genes,
                cdss,
                gene_lookup,
                args.intergenic_candidate_window,
                operon_map,
                args.tn5_structure_profile,
            )
            note = join_unique([chosen_note, annotation.note])
            main_rows.append(
                {
                    "Query_ID": record.query_id,
                    "Row_number": record.row_number,
                    "Gene name": annotation.gene_name,
                    "Gene code": annotation.gene_code,
                    "Insertion site": chosen.top_hit.insertion_site,
                    "Junction_boundary_0based": chosen.top_hit.junction_boundary_0based,
                    "Coordinate_definition": JUNCTION_COORDINATE_DEFINITION,
                    "Putative function or description": annotation.putative_function,
                    "Ontology_term": annotation.ontology_term,
                    "Go_tags": annotation.go_tags,
                    "Feature_relation": annotation.feature_relation,
                    "Candidate gene name": annotation.candidate_gene_name,
                    "Candidate gene code": annotation.candidate_gene_code,
                    "Candidate relation": annotation.candidate_relation,
                    "Candidate distance (bp)": annotation.candidate_distance_bp,
                    "Candidate strand": annotation.candidate_strand,
                    "Candidate priority": annotation.candidate_priority,
                    "Candidate function": annotation.candidate_function,
                    "Candidate ontology term": annotation.candidate_ontology_term,
                    "Candidate GO tags": annotation.candidate_go_tags,
                    "Candidate score": annotation.candidate_score,
                    "Candidate evidence": annotation.candidate_evidence,
                    "Candidate operon ID": annotation.candidate_operon_id,
                    "Candidate operon relation": annotation.candidate_operon_relation,
                    "Tn5 structure profile": annotation.tn5_structure_profile,
                    "Reference_seqid": chosen.top_hit.reference_seqid,
                    "Strand": chosen.top_hit.strand,
                    "Hit_status": chosen.status,
                    "Alignment_method": chosen.top_hit.alignment_method,
                    "Alignment_score": round(chosen.top_hit.alignment_score, 6),
                    "Score_gap": round(chosen.top_hit.score_gap, 6),
                    "Second_best_score": round(chosen.hits[1].alignment_score, 6) if len(chosen.hits) > 1 else "",
                    "Second_best_identity": round(chosen.hits[1].identity, 6) if len(chosen.hits) > 1 else "",
                    "Edit_distance": chosen.top_hit.edit_distance,
                    "Note": note,
                }
            )
            summary_counts[chosen.status] += 1
            if rescued:
                summary_counts["rescued_total"] += 1
            if annotation.feature_relation == "intergenic":
                summary_counts["intergenic"] += 1
                if annotation.candidate_gene_code or annotation.candidate_gene_name:
                    summary_counts["intergenic_with_nearby_candidate"] += 1
            else:
                summary_counts["gene_hits"] += 1
        else:
            main_rows.append(
                {
                    "Query_ID": record.query_id,
                    "Row_number": record.row_number,
                    "Gene name": "",
                    "Gene code": "",
                    "Insertion site": chosen.top_hit.insertion_site if chosen.top_hit else "",
                    "Junction_boundary_0based": chosen.top_hit.junction_boundary_0based if chosen.top_hit else "",
                    "Coordinate_definition": JUNCTION_COORDINATE_DEFINITION,
                    "Putative function or description": "",
                    "Ontology_term": "",
                    "Go_tags": "",
                    "Feature_relation": "",
                    "Reference_seqid": chosen.top_hit.reference_seqid if chosen.top_hit else "",
                    "Strand": chosen.top_hit.strand if chosen.top_hit else "",
                    "Hit_status": chosen.status,
                    "Alignment_method": chosen.top_hit.alignment_method if chosen.top_hit else aligner,
                    "Alignment_score": round(chosen.top_hit.alignment_score, 6) if chosen.top_hit else "",
                    "Score_gap": round(chosen.top_hit.score_gap, 6) if chosen.top_hit else "",
                    "Second_best_score": round(chosen.hits[1].alignment_score, 6) if len(chosen.hits) > 1 else "",
                    "Second_best_identity": round(chosen.hits[1].identity, 6) if len(chosen.hits) > 1 else "",
                    "Edit_distance": chosen.top_hit.edit_distance if chosen.top_hit else "",
                    "Note": chosen_note,
                }
            )
            summary_counts[chosen.status] += 1
            failed_rows.append(
                {
                    "Query_ID": record.query_id,
                    "Row_number": record.row_number,
                    "Category": chosen.status,
                    "Reason": chosen_note,
                }
            )

        raw_hits = chosen.hits or []
        if not raw_hits and chosen.top_hit is None:
            raw_match_rows.append(
                {
                    "Query_ID": record.query_id,
                    "Row_number": record.row_number,
                    "Candidate_source": chosen.candidate.source,
                    "Reference_seqid": "",
                    "Strand": "",
                    "Alignment_start": "",
                    "Alignment_end": "",
                    "Query_pos1_reference_coord": "",
                    "Insertion_site": "",
                    "Junction_boundary_0based": "",
                    "Identity": "",
                    "Query_coverage": "",
                    "Alignment_score": "",
                    "Score_gap": "",
                    "Edit_distance": "",
                    "Alignment_method": aligner,
                    "Alignment_length": len(chosen.candidate.sequence),
                    "Hit_count": 0,
                    "Hit_status": chosen.status,
                    "Note": chosen_note,
                }
            )
        else:
            for hit in raw_hits[:25]:
                raw_match_rows.append(
                    {
                        "Query_ID": hit.query_id,
                        "Row_number": hit.row_number,
                        "Candidate_source": hit.candidate_source,
                        "Reference_seqid": hit.reference_seqid,
                        "Strand": hit.strand,
                        "Alignment_start": hit.alignment_start,
                        "Alignment_end": hit.alignment_end,
                        "Query_pos1_reference_coord": hit.query_pos1_reference_coord,
                        "Insertion_site": hit.insertion_site,
                        "Junction_boundary_0based": hit.junction_boundary_0based,
                        "Identity": round(hit.identity, 6),
                        "Query_coverage": round(hit.query_coverage, 6),
                        "Alignment_score": round(hit.alignment_score, 6),
                        "Score_gap": round(hit.score_gap, 6),
                        "Edit_distance": hit.edit_distance,
                        "Alignment_method": hit.alignment_method,
                        "Alignment_length": hit.alignment_length,
                        "Hit_count": hit.hit_count,
                        "Hit_status": hit.hit_status,
                        "Note": hit.note,
                    }
                )

    main_headers = [
        "Query_ID",
        "Row_number",
        "Gene name",
        "Gene code",
        "Insertion site",
        "Junction_boundary_0based",
        "Coordinate_definition",
        "Putative function or description",
        "Ontology_term",
        "Go_tags",
        "Feature_relation",
        "Candidate gene name",
        "Candidate gene code",
        "Candidate relation",
        "Candidate distance (bp)",
        "Candidate strand",
        "Candidate priority",
        "Candidate function",
        "Candidate ontology term",
        "Candidate GO tags",
        "Candidate score",
        "Candidate evidence",
        "Candidate operon ID",
        "Candidate operon relation",
        "Tn5 structure profile",
        "Reference_seqid",
        "Strand",
        "Hit_status",
        "Alignment_method",
        "Alignment_score",
        "Score_gap",
        "Second_best_score",
        "Second_best_identity",
        "Edit_distance",
        "Note",
    ]
    raw_headers = [
        "Query_ID",
        "Row_number",
        "Candidate_source",
        "Reference_seqid",
        "Strand",
        "Alignment_start",
        "Alignment_end",
        "Query_pos1_reference_coord",
        "Insertion_site",
        "Junction_boundary_0based",
        "Identity",
        "Query_coverage",
        "Alignment_length",
        "Alignment_score",
        "Score_gap",
        "Edit_distance",
        "Alignment_method",
        "Hit_count",
        "Hit_status",
        "Note",
    ]
    failed_headers = ["Query_ID", "Row_number", "Category", "Reason"]
    try:
        from tn5cope.functional_clusters import (
            ASSIGNMENT_HEADERS as EXPANDED_ASSIGNMENT_HEADERS,
            EVENT_HEADERS as EXPANDED_EVENT_HEADERS,
            EXCLUSION_HEADERS as EXPANDED_EXCLUSION_HEADERS,
            TERM_SUMMARY_HEADERS as EXPANDED_TERM_SUMMARY_HEADERS,
            run_expanded_functional_cluster_analysis,
        )
        from tn5cope.position_clusters import (
            ASSIGNMENT_HEADERS as POSITION_ASSIGNMENT_HEADERS,
            SUMMARY_HEADERS as POSITION_SUMMARY_HEADERS,
            run_position_cluster_analysis,
        )
    except ModuleNotFoundError:  # pragma: no cover - supports direct script execution.
        from functional_clusters import (
            ASSIGNMENT_HEADERS as EXPANDED_ASSIGNMENT_HEADERS,
            EVENT_HEADERS as EXPANDED_EVENT_HEADERS,
            EXCLUSION_HEADERS as EXPANDED_EXCLUSION_HEADERS,
            TERM_SUMMARY_HEADERS as EXPANDED_TERM_SUMMARY_HEADERS,
            run_expanded_functional_cluster_analysis,
        )
        from position_clusters import (
            ASSIGNMENT_HEADERS as POSITION_ASSIGNMENT_HEADERS,
            SUMMARY_HEADERS as POSITION_SUMMARY_HEADERS,
            run_position_cluster_analysis,
        )

    expanded_summary = run_expanded_functional_cluster_analysis(
        main_rows,
        output_prefix,
    )
    position_summary = run_position_cluster_analysis(
        main_rows,
        expanded_summary["event_rows"],
        reference_seqid,
        len(genome),
        output_prefix,
        genes.get(reference_seqid, []),
    )
    summary_rows = [
        ["Metric", "Value"],
        ["total_queries", summary_counts["total"]],
        ["valid_sequences", summary_counts["valid_sequences"]],
        ["unique_hit", summary_counts["unique_hit"]],
        ["rescued_unique_hit", summary_counts["rescued_unique_hit"]],
        ["rescued_total", summary_counts["rescued_total"]],
        ["multiple_hits", summary_counts["multiple_hits"]],
        ["low_confidence", summary_counts["low_confidence"]],
        ["no_hit", summary_counts["no_hit"]],
        ["no_sequence", summary_counts["no_sequence"]],
        ["gene_hits", summary_counts["gene_hits"]],
        ["intergenic", summary_counts["intergenic"]],
        ["intergenic_with_nearby_candidate", summary_counts["intergenic_with_nearby_candidate"]],
        ["intergenic_candidate_window_bp", args.intergenic_candidate_window],
        ["operon_annotation", str(operon_path) if operon_path else "not_provided"],
        ["tn5_structure_profile", args.tn5_structure_profile],
        ["analysis_scope", "Preliminary candidate-gene analysis; experimental validation is required for causal claims."],
        ["expanded_direct_event_count", expanded_summary["metrics"]["direct_event_count"]],
        ["expanded_candidate_event_count", expanded_summary["metrics"]["candidate_event_count"]],
        ["expanded_direct_gene_count", expanded_summary["metrics"]["direct_gene_count"]],
        ["expanded_candidate_gene_count", expanded_summary["metrics"]["candidate_gene_count"]],
        ["expanded_overlapping_gene_count", expanded_summary["metrics"]["overlapping_gene_count"]],
        ["expanded_unique_gene_count", expanded_summary["metrics"]["unique_gene_count"]],
        ["expanded_named_cluster_count", expanded_summary["metrics"]["named_cluster_count"]],
        ["expanded_no_feature_gene_count", expanded_summary["metrics"]["no_feature_gene_count"]],
        ["expanded_candidate_exclusion_count", expanded_summary["metrics"]["candidate_exclusion_count"]],
        [
            "expanded_clustering_method_note",
            (
                "Reliable direct genes plus the uniquely highest-scoring high/moderate "
                "intergenic candidate; one equal-weight unit per unique gene; binary "
                "ontology/GO/product features, Jaccard distance, average linkage."
            ),
        ],
        ["unique_insertion_loci", position_summary["unique_insertion_loci"]],
        ["position_cluster_count", position_summary["position_cluster_count"]],
        ["position_noise_loci", position_summary["position_noise_loci"]],
        [
            "expanded_unique_insertion_loci",
            position_summary["expanded_unique_insertion_loci"],
        ],
        [
            "expanded_position_cluster_count",
            position_summary["expanded_position_cluster_count"],
        ],
        [
            "expanded_position_noise_loci",
            position_summary["expanded_position_noise_loci"],
        ],
        ["method_note", f"Alignment engine={aligner}; all informative seed sizes and positions are searched in the internal fallback, and standard local alignment is used when available."],
        ["rescue_method_note", "Second-pass rescue uses insertion-anchored terminal subfragments; multiple qualifying windows must agree on one junction coordinate and strand."],
        ["coordinate_definition", JUNCTION_COORDINATE_DEFINITION],
        [
            "position_cluster_method_note",
            (
                "Two circular DBSCAN analyses are reported: all unique high-confidence "
                "insertion coordinates and the Expanded candidate subset. The circle "
                "uses Expanded candidate PC spans/points and all high-confidence loci "
                "for the 50-kb density ring."
            ),
        ],
        ["identity_threshold_long", HIGH_IDENTITY_LONG],
        ["identity_threshold_short", HIGH_IDENTITY_SHORT],
        ["rescue_identity_threshold_long", RESCUE_IDENTITY_LONG],
        ["rescue_identity_threshold_short", RESCUE_IDENTITY_SHORT],
        ["min_query_coverage", MIN_QUERY_COVERAGE],
        ["top_hit_margin", TOP_HIT_MARGIN],
        ["top_hit_score_gap", TOP_HIT_SCORE_GAP],
    ]
    for reason, count in expanded_summary["metrics"]["candidate_exclusion_counts"].items():
        summary_rows.append([f"expanded_exclusion_{reason}", count])

    write_tsv(output_prefix.with_name(output_prefix.name + "_failed_or_ambiguous_hits.tsv"), failed_headers, failed_rows)
    write_tsv(output_prefix.with_name(output_prefix.name + "_main_results.tsv"), main_headers, main_rows)
    write_tsv(
        output_prefix.with_name(output_prefix.name + "_raw_matches.tsv"),
        raw_headers,
        raw_match_rows,
    )
    write_xlsx(
        output_prefix.with_suffix(".xlsx"),
        [
            ("main_results", rows_from_dicts(main_headers, main_rows)),
            ("raw_matches", rows_from_dicts(raw_headers, raw_match_rows)),
            ("failed_or_ambiguous", rows_from_dicts(failed_headers, failed_rows)),
            (
                "expanded_gene_events",
                rows_from_dicts(
                    EXPANDED_EVENT_HEADERS,
                    expanded_summary["event_rows"],
                ),
            ),
            (
                "expanded_clusters",
                rows_from_dicts(
                    EXPANDED_ASSIGNMENT_HEADERS,
                    expanded_summary["assignment_rows"],
                ),
            ),
            (
                "expanded_terms",
                rows_from_dicts(
                    EXPANDED_TERM_SUMMARY_HEADERS,
                    expanded_summary["term_rows"],
                ),
            ),
            (
                "candidate_exclusions",
                rows_from_dicts(
                    EXPANDED_EXCLUSION_HEADERS,
                    expanded_summary["exclusion_rows"],
                ),
            ),
            (
                "all_position_sites",
                rows_from_dicts(
                    POSITION_ASSIGNMENT_HEADERS,
                    position_summary["all_assignment_rows"],
                ),
            ),
            (
                "all_position_clusters",
                rows_from_dicts(
                    POSITION_SUMMARY_HEADERS,
                    position_summary["all_summary_rows"],
                ),
            ),
            (
                "expanded_position_sites",
                rows_from_dicts(
                    POSITION_ASSIGNMENT_HEADERS,
                    position_summary["expanded_assignment_rows"],
                ),
            ),
            (
                "expanded_position_PC",
                rows_from_dicts(
                    POSITION_SUMMARY_HEADERS,
                    position_summary["expanded_summary_rows"],
                ),
            ),
            ("summary", summary_rows),
        ],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
