import shutil

import pytest

from tn5cope.pipeline import (
    CandidateEvaluation,
    FragmentCandidate,
    Hit,
    SequenceRecord,
    build_kmer_index,
    choose_consistent_rescue_evaluation,
    junction_boundary_0based,
    map_fragment,
    map_fragment_external,
)


def test_all_informative_seed_sizes_are_searched() -> None:
    query = "ACGTGCAATG" * 8
    genome = "T" * 40 + query + "G" * 40 + query + "C" * 40
    record = SequenceRecord("Q001", 1, query, query)
    candidate = FragmentCandidate("full_query", query, "left", "test")
    indices = {k: build_kmer_index(genome, k) for k in (11, 9)}

    hits = map_fragment(record, candidate, "ctg", genome, indices, (11, 9))
    starts = {hit.alignment_start for hit in hits if hit.strand == "+"}
    assert len(starts) >= 2


def _rescue_evaluation(source: str, boundary: int, score: float) -> CandidateEvaluation:
    candidate = FragmentCandidate(source, "A" * 60, "left", "rescue")
    hit = Hit(
        query_id="Q001",
        row_number=1,
        candidate_source=source,
        reference_seqid="ctg",
        strand="+",
        alignment_start=boundary + 1,
        alignment_end=boundary + 60,
        query_pos1_reference_coord=boundary + 1,
        insertion_site=boundary + 1,
        identity=0.95,
        query_coverage=1.0,
        alignment_length=60,
        matched_bases=57,
        comparable_bases=60,
        alignment_score=score,
        junction_boundary_0based=boundary,
    )
    return CandidateEvaluation(candidate, [hit], "rescued_unique_hit", hit, "")


def test_conflicting_rescue_windows_are_ambiguous() -> None:
    chosen = choose_consistent_rescue_evaluation(
        [_rescue_evaluation("window150", 100, 100), _rescue_evaluation("window120", 300, 99)]
    )
    assert chosen is not None
    assert chosen.status == "multiple_hits"
    assert "conflicting junction" in chosen.note


def test_junction_coordinate_definition_is_explicit() -> None:
    assert junction_boundary_0based("+", "right", 101, 200) == 200
    assert junction_boundary_0based("-", "right", 101, 200) == 100


@pytest.mark.parametrize(
    ("aligner", "executable", "expected_method"),
    [
        ("minimap2", "minimap2", "minimap2"),
        ("blastn", "blastn", "blastn-short"),
        ("bwa", "bwa", "bwa"),
    ],
)
def test_installed_external_aligner_path(
    tmp_path, aligner: str, executable: str, expected_method: str
) -> None:
    if shutil.which(executable) is None:
        pytest.skip(f"{executable} is not installed")

    query = (
        "ACGTTGCAAGTCCGATGCTAACGGTTCAGATCGTACGACT"
        "GGAACCTTAGCGTTCGATACCGTTAAGC"
    )
    reference = "A" * 80 + query + "C" * 80
    genome_path = tmp_path / "reference.fna"
    genome_path.write_text(f">ctg\n{reference}\n", encoding="ascii")
    record = SequenceRecord("Q_external", 1, query, query)
    candidate = FragmentCandidate("full_query", query, "left", "integration test")

    hits = map_fragment_external(record, candidate, "ctg", genome_path, aligner)

    assert hits
    assert hits[0].alignment_start == 81
    assert hits[0].identity == pytest.approx(1.0)
    assert hits[0].query_coverage == pytest.approx(1.0)
    assert hits[0].alignment_method == expected_method
