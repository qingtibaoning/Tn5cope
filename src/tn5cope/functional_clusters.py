from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist

try:
    from tn5cope.position_clusters import setup_times_new_roman
except ModuleNotFoundError:  # pragma: no cover - supports direct script execution.
    from position_clusters import setup_times_new_roman


VALID_HIT_STATUSES = {"unique_hit", "rescued_unique_hit"}
DIRECT_FEATURE_RELATIONS = {"CDS", "gene_non_CDS"}
ELIGIBLE_CANDIDATE_PRIORITIES = {"high", "moderate"}
CANDIDATE_SEPARATOR = " | "
CANDIDATE_ANNOTATION_SEPARATOR = " || "
FUNCTIONAL_MIN_TERM_FREQUENCY = 2
FUNCTIONAL_MAX_TERM_PREVALENCE = 0.70
FUNCTIONAL_CLUSTER_CUT_DISTANCE = 0.70
FUNCTIONAL_MIN_CLUSTER_SIZE = 2
PLOT_DPI = 220

PRODUCT_STOPWORDS = {
    "protein",
    "putative",
    "family",
    "domain",
    "containing",
    "hypothetical",
    "probable",
    "predicted",
    "conserved",
    "subunit",
    "chain",
    "component",
    "dependent",
    "enzyme",
    "like",
    "related",
    "associated",
    "system",
    "type",
}

EVENT_HEADERS = [
    "Query_ID",
    "Row_number",
    "Insertion site",
    "Reference_seqid",
    "Hit_status",
    "Feature_relation",
    "Selected gene name",
    "Selected gene code",
    "Selected gene ID",
    "Gene source",
    "Evidence tier",
    "Candidate relation",
    "Candidate distance (bp)",
    "Candidate strand",
    "Candidate score",
    "Candidate evidence",
    "Candidate operon ID",
    "Candidate operon relation",
    "Product",
    "Ontology_term",
    "Go_tags",
    "Selection rule",
]

EXCLUSION_HEADERS = [
    "Query_ID",
    "Row_number",
    "Insertion site",
    "Reference_seqid",
    "Hit_status",
    "Feature_relation",
    "Exclusion reason",
    "Candidate count",
    "Top candidate score",
    "Candidate gene names",
    "Candidate gene codes",
    "Candidate priorities",
    "Candidate scores",
    "Details",
]

ASSIGNMENT_HEADERS = [
    "cluster_id",
    "cluster_size",
    "gene_id",
    "gene_name",
    "gene_code",
    "gene_sources",
    "evidence_tiers",
    "source_label",
    "clustering_weight",
    "supporting_events",
    "direct_supporting_events",
    "candidate_supporting_events",
    "query_ids",
    "insertion_sites",
    "feature_relations",
    "retained_feature_count",
    "retained_features",
    "products",
    "ontology_terms",
    "go_names",
    "hit_statuses",
]

TERM_SUMMARY_HEADERS = [
    "cluster_id",
    "cluster_size",
    "term_rank",
    "feature_type",
    "term",
    "gene_count",
    "cluster_prevalence",
    "member_genes",
]

SOURCE_ORDER = {
    "direct": 0,
    "candidate_high": 1,
    "candidate_moderate": 2,
}
EVIDENCE_ORDER = {"primary": 0, "high": 1, "moderate": 2}


@dataclass
class ExpandedGeneRecord:
    gene_id: str
    gene_names: set[str] = field(default_factory=set)
    gene_codes: set[str] = field(default_factory=set)
    events: List[Dict[str, object]] = field(default_factory=list)
    products: set[str] = field(default_factory=set)
    ontology_terms: set[str] = field(default_factory=set)
    go_names: set[str] = field(default_factory=set)
    hit_statuses: set[str] = field(default_factory=set)
    feature_relations: set[str] = field(default_factory=set)
    gene_sources: set[str] = field(default_factory=set)
    evidence_tiers: set[str] = field(default_factory=set)

    @property
    def gene_name(self) -> str:
        return sorted(self.gene_names)[0] if self.gene_names else ""

    @property
    def gene_code(self) -> str:
        return sorted(self.gene_codes)[0] if self.gene_codes else ""


def output_path(output_prefix: Path, suffix: str) -> Path:
    return output_prefix.with_name(output_prefix.name + suffix)


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _split_semicolon(value: object) -> List[str]:
    return [part.strip() for part in _text(value).split(";") if part.strip()]


def parse_go_names(go_tags: object) -> List[str]:
    names: List[str] = []
    for item in _split_semicolon(go_tags):
        name = item.split("|", 1)[0].strip().lower()
        if name:
            names.append(name)
    return names


def tokenize_product(text: str) -> List[str]:
    clean = text.lower().replace("/", " ").replace("-", " ")
    tokens = re.findall(r"[a-z0-9]+", clean)
    return [
        token
        for token in tokens
        if len(token) >= 3 and token not in PRODUCT_STOPWORDS and not token.isdigit()
    ]


def normalized_product_phrase(text: str) -> str:
    return " ".join(tokenize_product(text))


def _candidate_values(
    row: Dict[str, object],
    fields: Sequence[str],
    separator: str = CANDIDATE_SEPARATOR,
    expected_count: int | None = None,
) -> Tuple[Dict[str, List[str]], int, str]:
    split_values: Dict[str, List[str]] = {}
    nonempty_lengths: List[int] = []
    for field_name in fields:
        raw_value = row.get(field_name)
        raw = "" if raw_value is None else str(raw_value)
        if not raw.strip():
            split_values[field_name] = []
            continue
        values = [part.strip() for part in raw.split(separator)]
        split_values[field_name] = values
        nonempty_lengths.append(len(values))

    candidate_count = expected_count if expected_count is not None else max(nonempty_lengths, default=0)
    if candidate_count == 0:
        return split_values, 0, ""

    for field_name in fields:
        values = split_values[field_name]
        if not values:
            split_values[field_name] = [""] * candidate_count
        elif len(values) != candidate_count:
            return (
                split_values,
                candidate_count,
                f"{field_name} has {len(values)} values; expected {candidate_count}",
            )
    return split_values, candidate_count, ""


def _base_exclusion(row: Dict[str, object], reason: str, details: str = "") -> Dict[str, object]:
    return {
        "Query_ID": row.get("Query_ID", ""),
        "Row_number": row.get("Row_number", ""),
        "Insertion site": row.get("Insertion site", ""),
        "Reference_seqid": row.get("Reference_seqid", ""),
        "Hit_status": row.get("Hit_status", ""),
        "Feature_relation": row.get("Feature_relation", ""),
        "Exclusion reason": reason,
        "Candidate count": "",
        "Top candidate score": "",
        "Candidate gene names": row.get("Candidate gene name", ""),
        "Candidate gene codes": row.get("Candidate gene code", ""),
        "Candidate priorities": row.get("Candidate priority", ""),
        "Candidate scores": row.get("Candidate score", ""),
        "Details": details,
    }


def _direct_event(row: Dict[str, object]) -> Tuple[Dict[str, object] | None, Dict[str, object] | None]:
    gene_name = _text(row.get("Gene name"))
    gene_code = _text(row.get("Gene code"))
    gene_id = gene_code or gene_name
    if not gene_id:
        return None, _base_exclusion(
            row,
            "direct_gene_missing_identifier",
            "A direct gene interval was reported without a gene code or gene name.",
        )

    product = _text(row.get("Putative function or description")) or _text(
        row.get("Candidate function")
    )
    ontology = _text(row.get("Ontology_term")) or _text(
        row.get("Candidate ontology term")
    )
    go_tags = _text(row.get("Go_tags")) or _text(row.get("Candidate GO tags"))
    return (
        {
            "Query_ID": row.get("Query_ID", ""),
            "Row_number": row.get("Row_number", ""),
            "Insertion site": row.get("Insertion site", ""),
            "Reference_seqid": row.get("Reference_seqid", ""),
            "Hit_status": row.get("Hit_status", ""),
            "Feature_relation": row.get("Feature_relation", ""),
            "Selected gene name": gene_name,
            "Selected gene code": gene_code,
            "Selected gene ID": gene_id,
            "Gene source": "direct",
            "Evidence tier": "primary",
            "Candidate relation": row.get("Candidate relation", ""),
            "Candidate distance (bp)": 0,
            "Candidate strand": row.get("Candidate strand", ""),
            "Candidate score": row.get("Candidate score", ""),
            "Candidate evidence": row.get("Candidate evidence", ""),
            "Candidate operon ID": row.get("Candidate operon ID", ""),
            "Candidate operon relation": row.get("Candidate operon relation", ""),
            "Product": product,
            "Ontology_term": ontology,
            "Go_tags": go_tags,
            "Selection rule": "reliable_direct_gene_hit",
        },
        None,
    )


def _candidate_event(
    row: Dict[str, object],
) -> Tuple[Dict[str, object] | None, Dict[str, object] | None]:
    core_fields = [
        "Candidate gene name",
        "Candidate gene code",
        "Candidate relation",
        "Candidate distance (bp)",
        "Candidate strand",
        "Candidate priority",
        "Candidate function",
        "Candidate score",
        "Candidate evidence",
        "Candidate operon ID",
        "Candidate operon relation",
    ]
    values, candidate_count, error = _candidate_values(row, core_fields)
    if candidate_count == 0:
        return None, _base_exclusion(
            row,
            "no_candidates",
            "No gene was reported within the maximum intergenic candidate window.",
        )
    if error:
        exclusion = _base_exclusion(row, "candidate_columns_misaligned", error)
        exclusion["Candidate count"] = candidate_count
        return None, exclusion

    annotation_fields = ["Candidate ontology term", "Candidate GO tags"]
    annotations, _, annotation_error = _candidate_values(
        row,
        annotation_fields,
        separator=CANDIDATE_ANNOTATION_SEPARATOR,
        expected_count=candidate_count,
    )
    if annotation_error:
        exclusion = _base_exclusion(
            row, "candidate_columns_misaligned", annotation_error
        )
        exclusion["Candidate count"] = candidate_count
        return None, exclusion

    scores: List[float] = []
    try:
        for raw_score in values["Candidate score"]:
            scores.append(float(raw_score))
    except ValueError:
        exclusion = _base_exclusion(
            row,
            "invalid_candidate_score",
            "Every reported candidate must have a numeric score.",
        )
        exclusion["Candidate count"] = candidate_count
        return None, exclusion

    top_score = max(scores)
    top_indices = [index for index, score in enumerate(scores) if score == top_score]
    if len(top_indices) != 1:
        exclusion = _base_exclusion(
            row,
            "top_score_tied",
            "More than one candidate shares the highest score.",
        )
        exclusion["Candidate count"] = candidate_count
        exclusion["Top candidate score"] = top_score
        return None, exclusion

    selected_index = top_indices[0]
    priority = values["Candidate priority"][selected_index].lower()
    if priority not in ELIGIBLE_CANDIDATE_PRIORITIES:
        exclusion = _base_exclusion(
            row,
            "top_priority_not_eligible",
            "Only uniquely top-ranked high or moderate intergenic candidates are included.",
        )
        exclusion["Candidate count"] = candidate_count
        exclusion["Top candidate score"] = top_score
        return None, exclusion

    gene_name = values["Candidate gene name"][selected_index]
    gene_code = values["Candidate gene code"][selected_index]
    gene_id = gene_code or gene_name
    if not gene_id:
        exclusion = _base_exclusion(
            row,
            "top_candidate_missing_identifier",
            "The uniquely top-ranked candidate has neither a gene code nor a gene name.",
        )
        exclusion["Candidate count"] = candidate_count
        exclusion["Top candidate score"] = top_score
        return None, exclusion

    return (
        {
            "Query_ID": row.get("Query_ID", ""),
            "Row_number": row.get("Row_number", ""),
            "Insertion site": row.get("Insertion site", ""),
            "Reference_seqid": row.get("Reference_seqid", ""),
            "Hit_status": row.get("Hit_status", ""),
            "Feature_relation": row.get("Feature_relation", ""),
            "Selected gene name": gene_name,
            "Selected gene code": gene_code,
            "Selected gene ID": gene_id,
            "Gene source": f"candidate_{priority}",
            "Evidence tier": priority,
            "Candidate relation": values["Candidate relation"][selected_index],
            "Candidate distance (bp)": values["Candidate distance (bp)"][selected_index],
            "Candidate strand": values["Candidate strand"][selected_index],
            "Candidate score": values["Candidate score"][selected_index],
            "Candidate evidence": values["Candidate evidence"][selected_index],
            "Candidate operon ID": values["Candidate operon ID"][selected_index],
            "Candidate operon relation": values["Candidate operon relation"][selected_index],
            "Product": values["Candidate function"][selected_index],
            "Ontology_term": annotations["Candidate ontology term"][selected_index],
            "Go_tags": annotations["Candidate GO tags"][selected_index],
            "Selection rule": "unique_top_score_and_high_or_moderate_priority",
        },
        None,
    )


def select_expanded_gene_events(
    rows: Sequence[Dict[str, object]],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    events: List[Dict[str, object]] = []
    exclusions: List[Dict[str, object]] = []
    seen_query_gene: set[Tuple[str, str]] = set()

    for row in rows:
        if _text(row.get("Hit_status")) not in VALID_HIT_STATUSES:
            continue
        relation = _text(row.get("Feature_relation"))
        if relation in DIRECT_FEATURE_RELATIONS:
            event, exclusion = _direct_event(row)
        elif relation == "intergenic":
            event, exclusion = _candidate_event(row)
        else:
            event = None
            exclusion = _base_exclusion(
                row,
                "unsupported_feature_relation",
                f"Reliable hit has unsupported feature relation {relation!r}.",
            )

        if exclusion is not None:
            exclusions.append(exclusion)
        if event is None:
            continue
        key = (_text(event["Query_ID"]), _text(event["Selected gene ID"]))
        if key in seen_query_gene:
            continue
        seen_query_gene.add(key)
        events.append(event)

    events.sort(key=lambda item: (_text(item["Query_ID"]), _text(item["Selected gene ID"])))
    exclusions.sort(key=lambda item: (_text(item["Query_ID"]), _text(item["Exclusion reason"])))
    return events, exclusions


def _sorted_join(values: Iterable[str], order: Dict[str, int] | None = None) -> str:
    unique = {value for value in values if value}
    if order is None:
        ordered = sorted(unique)
    else:
        ordered = sorted(unique, key=lambda value: (order.get(value, 999), value))
    return "; ".join(ordered)


def aggregate_gene_records(events: Sequence[Dict[str, object]]) -> List[ExpandedGeneRecord]:
    grouped: Dict[str, ExpandedGeneRecord] = {}
    for event in events:
        gene_id = _text(event.get("Selected gene ID"))
        record = grouped.setdefault(gene_id, ExpandedGeneRecord(gene_id=gene_id))
        gene_name = _text(event.get("Selected gene name"))
        gene_code = _text(event.get("Selected gene code"))
        if gene_name:
            record.gene_names.add(gene_name)
        if gene_code:
            record.gene_codes.add(gene_code)
        record.events.append(dict(event))
        product = _text(event.get("Product"))
        if product:
            record.products.add(product)
        record.ontology_terms.update(
            term.lower() for term in _split_semicolon(event.get("Ontology_term"))
        )
        record.go_names.update(parse_go_names(event.get("Go_tags")))
        record.hit_statuses.add(_text(event.get("Hit_status")))
        record.feature_relations.add(_text(event.get("Feature_relation")))
        record.gene_sources.add(_text(event.get("Gene source")))
        record.evidence_tiers.add(_text(event.get("Evidence tier")))
    return [grouped[gene_id] for gene_id in sorted(grouped)]


def build_functional_feature_map(
    records: Sequence[ExpandedGeneRecord],
) -> Tuple[Dict[str, List[str]], Dict[str, int]]:
    gene_features: Dict[str, List[str]] = {}
    feature_counts: Counter[str] = Counter()
    n_genes = len(records)
    max_prevalence_count = max(
        1, math.floor(FUNCTIONAL_MAX_TERM_PREVALENCE * n_genes)
    )

    for record in records:
        features: set[str] = set()
        features.update(f"ontology:{term}" for term in record.ontology_terms if term)
        features.update(f"go:{name}" for name in record.go_names if name)
        for product in record.products:
            phrase = normalized_product_phrase(product)
            if phrase:
                features.add(f"product_phrase:{phrase}")
            features.update(
                f"product_token:{token}" for token in tokenize_product(product)
            )
        gene_features[record.gene_id] = sorted(features)
        feature_counts.update(features)

    retained = {
        feature
        for feature, count in feature_counts.items()
        if count >= FUNCTIONAL_MIN_TERM_FREQUENCY
        and count <= max_prevalence_count
    }
    filtered = {
        gene_id: [feature for feature in features if feature in retained]
        for gene_id, features in gene_features.items()
    }
    return filtered, {feature: feature_counts[feature] for feature in retained}


def _source_label(sources: Iterable[str]) -> str:
    labels = {
        "direct": "D",
        "candidate_high": "H",
        "candidate_moderate": "M",
    }
    return "+".join(
        labels[source]
        for source in sorted(
            {source for source in sources if source},
            key=lambda source: (SOURCE_ORDER.get(source, 999), source),
        )
        if source in labels
    )


def functional_cluster_analysis(records: Sequence[ExpandedGeneRecord]) -> Dict[str, object]:
    gene_features, retained_counts = build_functional_feature_map(records)
    clustered_records = [
        record for record in records if gene_features.get(record.gene_id)
    ]
    unclustered_records = [
        record for record in records if not gene_features.get(record.gene_id)
    ]
    ordered_features = sorted(
        retained_counts, key=lambda feature: (-retained_counts[feature], feature)
    )
    gene_order = sorted(record.gene_id for record in clustered_records)
    matrix = np.zeros((len(gene_order), len(ordered_features)), dtype=int)
    for row_index, gene_id in enumerate(gene_order):
        feature_set = set(gene_features[gene_id])
        for column_index, feature in enumerate(ordered_features):
            matrix[row_index, column_index] = int(feature in feature_set)

    if len(gene_order) >= 2:
        distance_vector = pdist(matrix, metric="jaccard")
        linkage_matrix = linkage(distance_vector, method="average")
        raw_labels = fcluster(
            linkage_matrix,
            t=FUNCTIONAL_CLUSTER_CUT_DISTANCE,
            criterion="distance",
        )
        leaf_order = list(dendrogram(linkage_matrix, no_plot=True)["leaves"])
    else:
        linkage_matrix = None
        raw_labels = np.array([1] * len(gene_order), dtype=int)
        leaf_order = list(range(len(gene_order)))

    cluster_members: Dict[int, List[str]] = defaultdict(list)
    for gene_id, raw_label in zip(gene_order, raw_labels):
        cluster_members[int(raw_label)].append(gene_id)

    named_groups = sorted(
        (
            (raw_label, sorted(members))
            for raw_label, members in cluster_members.items()
            if len(members) >= FUNCTIONAL_MIN_CLUSTER_SIZE
        ),
        key=lambda item: (-len(item[1]), item[1]),
    )
    named_cluster_ids = {
        raw_label: f"FC{index:02d}"
        for index, (raw_label, _) in enumerate(named_groups, start=1)
    }
    record_by_id = {record.gene_id: record for record in records}
    assignments: List[Dict[str, object]] = []

    for gene_id, raw_label in zip(gene_order, raw_labels):
        record = record_by_id[gene_id]
        members = cluster_members[int(raw_label)]
        cluster_id = named_cluster_ids.get(int(raw_label), "SINGLETON")
        assignments.append(
            _assignment_row(
                record,
                cluster_id,
                len(members),
                gene_features[gene_id],
            )
        )
    for record in unclustered_records:
        assignments.append(
            _assignment_row(record, "NO_RETAINED_FEATURES", 1, [])
        )

    term_summary: List[Dict[str, object]] = []
    for raw_label, members in named_groups:
        cluster_id = named_cluster_ids[raw_label]
        term_counter: Counter[str] = Counter()
        for gene_id in members:
            term_counter.update(gene_features[gene_id])
        for rank, (feature, count) in enumerate(
            sorted(term_counter.items(), key=lambda item: (-item[1], item[0])),
            start=1,
        ):
            feature_type, term = feature.split(":", 1)
            term_summary.append(
                {
                    "cluster_id": cluster_id,
                    "cluster_size": len(members),
                    "term_rank": rank,
                    "feature_type": feature_type,
                    "term": term,
                    "gene_count": count,
                    "cluster_prevalence": round(count / len(members), 6),
                    "member_genes": "; ".join(members),
                }
            )

    cluster_sort = {"NO_RETAINED_FEATURES": 998, "SINGLETON": 999}
    assignments.sort(
        key=lambda row: (
            cluster_sort.get(_text(row["cluster_id"]), 0),
            _text(row["cluster_id"]),
            _text(row["gene_id"]),
        )
    )
    return {
        "matrix": matrix,
        "feature_names": ordered_features,
        "gene_order": gene_order,
        "leaf_order": leaf_order,
        "linkage_matrix": linkage_matrix,
        "assignments": assignments,
        "term_summary": term_summary,
        "retained_counts": retained_counts,
        "named_cluster_count": len(named_groups),
    }


def _assignment_row(
    record: ExpandedGeneRecord,
    cluster_id: str,
    cluster_size: int,
    retained_features: Sequence[str],
) -> Dict[str, object]:
    query_ids = sorted(
        {_text(event.get("Query_ID")) for event in record.events if event.get("Query_ID")}
    )
    insertion_sites = sorted(
        {
            int(event["Insertion site"])
            for event in record.events
            if event.get("Insertion site") not in ("", None)
        }
    )
    direct_count = sum(
        1 for event in record.events if event.get("Gene source") == "direct"
    )
    candidate_count = len(record.events) - direct_count
    return {
        "cluster_id": cluster_id,
        "cluster_size": cluster_size,
        "gene_id": record.gene_id,
        "gene_name": record.gene_name,
        "gene_code": record.gene_code,
        "gene_sources": _sorted_join(record.gene_sources, SOURCE_ORDER),
        "evidence_tiers": _sorted_join(record.evidence_tiers, EVIDENCE_ORDER),
        "source_label": _source_label(record.gene_sources),
        "clustering_weight": 1,
        "supporting_events": len(record.events),
        "direct_supporting_events": direct_count,
        "candidate_supporting_events": candidate_count,
        "query_ids": "; ".join(query_ids),
        "insertion_sites": "; ".join(str(site) for site in insertion_sites),
        "feature_relations": _sorted_join(record.feature_relations),
        "retained_feature_count": len(retained_features),
        "retained_features": "; ".join(retained_features),
        "products": _sorted_join(record.products),
        "ontology_terms": _sorted_join(record.ontology_terms),
        "go_names": _sorted_join(record.go_names),
        "hit_statuses": _sorted_join(record.hit_statuses),
    }


def write_tsv(
    path: Path,
    headers: Sequence[str],
    rows: Sequence[Dict[str, object]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(headers), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _display_feature(feature: str) -> str:
    prefix, term = feature.split(":", 1)
    abbreviations = {
        "ontology": "ONT",
        "go": "GO",
        "product_phrase": "PHRASE",
        "product_token": "TOKEN",
    }
    text = f"{abbreviations.get(prefix, prefix.upper())}: {term}"
    return text if len(text) <= 48 else text[:45] + "..."


def plot_functional_heatmap(
    analysis: Dict[str, object],
    records: Sequence[ExpandedGeneRecord],
    png_path: Path,
    pdf_path: Path,
) -> None:
    setup_times_new_roman()
    gene_order = list(analysis["gene_order"])
    feature_names = list(analysis["feature_names"])
    matrix = np.asarray(analysis["matrix"])
    linkage_matrix = analysis["linkage_matrix"]
    record_by_id = {record.gene_id: record for record in records}

    if not gene_order or not feature_names:
        fig, axis = plt.subplots(figsize=(10, 4.8))
        axis.axis("off")
        axis.text(
            0.5,
            0.62,
            "Expanded candidate functional clustering",
            ha="center",
            va="center",
            fontsize=17,
            fontweight="bold",
        )
        axis.text(
            0.5,
            0.42,
            (
                f"{len(records)} genes were selected, but no functional term passed "
                f"the frequency ({FUNCTIONAL_MIN_TERM_FREQUENCY}) and prevalence "
                f"({FUNCTIONAL_MAX_TERM_PREVALENCE:.0%}) filters."
            ),
            ha="center",
            va="center",
            fontsize=11,
            color="#4b5563",
            wrap=True,
        )
        fig.savefig(png_path, dpi=PLOT_DPI, bbox_inches="tight", facecolor="white")
        fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return

    height = max(7.0, min(18.0, 0.28 * len(gene_order) + 2.5))
    width = max(12.0, min(24.0, 0.36 * len(feature_names) + 7.0))
    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=(1.8, max(3.5, 0.30 * len(feature_names))),
        wspace=0.03,
    )
    dendrogram_axis = fig.add_subplot(grid[0, 0])
    heatmap_axis = fig.add_subplot(grid[0, 1])

    if linkage_matrix is not None:
        dendrogram_data = dendrogram(
            linkage_matrix,
            orientation="left",
            no_labels=True,
            color_threshold=FUNCTIONAL_CLUSTER_CUT_DISTANCE,
            above_threshold_color="#7c8798",
            ax=dendrogram_axis,
        )
        leaf_order = list(dendrogram_data["leaves"])
        dendrogram_axis.axvline(
            FUNCTIONAL_CLUSTER_CUT_DISTANCE,
            color="#d1495b",
            linestyle="--",
            linewidth=1.2,
        )
    else:
        leaf_order = [0]
        dendrogram_axis.set_ylim(-0.5, 0.5)
    dendrogram_axis.set_xlabel("Jaccard distance", fontsize=9)
    dendrogram_axis.set_yticks([])
    dendrogram_axis.spines[["top", "right", "left"]].set_visible(False)

    ordered_matrix = matrix[leaf_order, :]
    ordered_genes = [gene_order[index] for index in leaf_order]
    heatmap_axis.imshow(
        ordered_matrix,
        aspect="auto",
        interpolation="nearest",
        origin="lower",
        cmap=matplotlib.colors.ListedColormap(["#f1f5f9", "#2563eb"]),
        vmin=0,
        vmax=1,
    )
    y_labels = [
        f"{gene_id} [{_source_label(record_by_id[gene_id].gene_sources)}]"
        for gene_id in ordered_genes
    ]
    heatmap_axis.set_yticks(np.arange(len(ordered_genes)), labels=y_labels, fontsize=8)
    heatmap_axis.yaxis.tick_right()
    heatmap_axis.set_xticks(
        np.arange(len(feature_names)),
        labels=[_display_feature(feature) for feature in feature_names],
        rotation=60,
        ha="right",
        fontsize=7,
    )
    heatmap_axis.tick_params(axis="both", length=0)
    heatmap_axis.set_xticks(
        np.arange(-0.5, len(feature_names), 1), minor=True
    )
    heatmap_axis.set_yticks(np.arange(-0.5, len(ordered_genes), 1), minor=True)
    heatmap_axis.grid(which="minor", color="white", linewidth=0.55)
    heatmap_axis.set_ylabel("Gene [evidence source]", fontsize=9)
    heatmap_axis.yaxis.set_label_position("right")
    for spine in heatmap_axis.spines.values():
        spine.set_visible(False)

    fig.suptitle(
        "Expanded candidate functional clustering",
        x=0.02,
        y=0.995,
        ha="left",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.02,
        0.967,
        (
            "Direct hits and uniquely top-ranked high/moderate intergenic candidates; "
            "equal gene weights. D = direct, H = high, M = moderate."
        ),
        ha="left",
        va="top",
        fontsize=9.5,
        color="#4b5563",
    )
    fig.subplots_adjust(top=0.92, bottom=0.25, left=0.04, right=0.91)
    fig.savefig(png_path, dpi=PLOT_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_methods(path: Path, metrics: Dict[str, object]) -> None:
    text = f"""Expanded candidate functional clustering

Scope
Only reliable unique_hit and rescued_unique_hit events are eligible. Direct CDS and
gene_non_CDS hits are included. An intergenic event is included only when exactly one
candidate has the highest candidate score and that candidate has high or moderate
priority. Tied top scores, low/context-only priorities, missing identifiers, malformed
candidate columns, and absent candidates are excluded with an explicit reason.

Unit and weighting
The clustering unit is one unique gene identifier (gene code preferred, otherwise gene
name). Repeated insertion events are aggregated. A gene supported by both direct and
candidate evidence occurs once, with all source labels retained. Every gene has
clustering_weight=1; direct and candidate-supported genes therefore have equal weight.

Features and clustering
Binary features are built from ontology terms, GO names, normalized product phrases,
and normalized product tokens. Product stopwords are removed. Features must occur in
at least {FUNCTIONAL_MIN_TERM_FREQUENCY} genes and no more than
{FUNCTIONAL_MAX_TERM_PREVALENCE:.0%} of selected genes. Gene distances are Jaccard
distances, hierarchical linkage is average linkage, and the dendrogram is cut at
distance {FUNCTIONAL_CLUSTER_CUT_DISTANCE:.2f}. Only groups with at least
{FUNCTIONAL_MIN_CLUSTER_SIZE} genes receive FC cluster identifiers. Singletons and
genes without retained features remain explicitly reported.

Interpretation
The 500-bp intergenic window is a maximum candidate reporting window, not a biological
effect threshold. Candidate ranking incorporates transcriptional direction, distance
tier, available operon information, and the selected Tn5 structure profile. These
clusters are preliminary hypothesis-generating results and do not establish causality.

Run summary
Direct events: {metrics['direct_event_count']}
Candidate events: {metrics['candidate_event_count']}
Unique expanded genes: {metrics['unique_gene_count']}
Genes supported by both direct and candidate events: {metrics['overlapping_gene_count']}
Named functional clusters: {metrics['named_cluster_count']}
Genes without retained functional features: {metrics['no_feature_gene_count']}
"""
    path.write_text(text, encoding="utf-8")


def run_expanded_functional_cluster_analysis(
    main_rows: Sequence[Dict[str, object]],
    output_prefix: Path,
) -> Dict[str, object]:
    event_rows, exclusion_rows = select_expanded_gene_events(main_rows)
    records = aggregate_gene_records(event_rows)
    analysis = functional_cluster_analysis(records)
    assignment_rows = list(analysis["assignments"])
    term_rows = list(analysis["term_summary"])

    direct_event_count = sum(
        1 for event in event_rows if event["Gene source"] == "direct"
    )
    candidate_event_count = len(event_rows) - direct_event_count
    direct_genes = {
        _text(event["Selected gene ID"])
        for event in event_rows
        if event["Gene source"] == "direct"
    }
    candidate_genes = {
        _text(event["Selected gene ID"])
        for event in event_rows
        if str(event["Gene source"]).startswith("candidate_")
    }
    exclusion_counts = Counter(
        _text(exclusion["Exclusion reason"]) for exclusion in exclusion_rows
    )
    metrics: Dict[str, object] = {
        "direct_event_count": direct_event_count,
        "candidate_event_count": candidate_event_count,
        "direct_gene_count": len(direct_genes),
        "candidate_gene_count": len(candidate_genes),
        "overlapping_gene_count": len(direct_genes & candidate_genes),
        "unique_gene_count": len(records),
        "named_cluster_count": analysis["named_cluster_count"],
        "no_feature_gene_count": sum(
            1
            for assignment in assignment_rows
            if assignment["cluster_id"] == "NO_RETAINED_FEATURES"
        ),
        "candidate_exclusion_count": len(exclusion_rows),
        "candidate_exclusion_counts": dict(sorted(exclusion_counts.items())),
    }

    write_tsv(
        output_path(output_prefix, "_expanded_gene_events.tsv"),
        EVENT_HEADERS,
        event_rows,
    )
    write_tsv(
        output_path(output_prefix, "_expanded_candidate_exclusions.tsv"),
        EXCLUSION_HEADERS,
        exclusion_rows,
    )
    write_tsv(
        output_path(output_prefix, "_expanded_functional_cluster_assignments.tsv"),
        ASSIGNMENT_HEADERS,
        assignment_rows,
    )
    write_tsv(
        output_path(output_prefix, "_expanded_functional_cluster_term_summary.tsv"),
        TERM_SUMMARY_HEADERS,
        term_rows,
    )
    plot_functional_heatmap(
        analysis,
        records,
        output_path(output_prefix, "_expanded_functional_clustering_heatmap.png"),
        output_path(output_prefix, "_expanded_functional_clustering_heatmap.pdf"),
    )
    write_methods(
        output_path(output_prefix, "_expanded_functional_clustering_methods.txt"),
        metrics,
    )
    return {
        "event_rows": event_rows,
        "exclusion_rows": exclusion_rows,
        "assignment_rows": assignment_rows,
        "term_rows": term_rows,
        "metrics": metrics,
    }
