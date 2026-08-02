from __future__ import annotations

from pathlib import Path

from tn5cope.pipeline import (
    CDSFeature,
    GeneFeature,
    annotate_insertion,
    build_gene_lookup,
    load_operon_map,
    outside_gene_relation,
)


def _gene(code: str, start: int, end: int, strand: str, name: str = "") -> GeneFeature:
    return GeneFeature(
        seqid="ctg",
        start=start,
        end=end,
        strand=strand,
        gene_name=name,
        gene_code=code,
        locus_tag=f"LT_{code}",
        gene_id=f"ID_{code}",
        attributes={},
    )


def _annotate(
    coord: int,
    genes_list: list[GeneFeature],
    *,
    operon_map: dict[str, str] | None = None,
    profile: str = "unknown",
):
    genes = {"ctg": genes_list}
    return annotate_insertion(
        "ctg",
        coord,
        genes,
        {"ctg": []},
        build_gene_lookup(genes),
        500,
        operon_map,
        profile,
    )


def test_transcriptional_side_is_strand_aware() -> None:
    plus = _gene("plus", 200, 300, "+")
    minus = _gene("minus", 200, 300, "-")

    assert outside_gene_relation(plus, 190)[:2] == ("5prime_upstream", 10)
    assert outside_gene_relation(minus, 190)[:2] == ("3prime_downstream", 10)
    assert outside_gene_relation(plus, 310)[:2] == ("3prime_downstream", 10)
    assert outside_gene_relation(minus, 310)[:2] == ("5prime_upstream", 10)


def test_intergenic_candidates_are_ranked_by_transcription_direction() -> None:
    gene_a = _gene("geneA", 100, 200, "-", "geneA")
    gene_b = _gene("geneB", 264, 400, "-", "geneB")
    result = _annotate(252, [gene_a, gene_b])

    assert result.feature_relation == "intergenic"
    assert result.candidate_gene_code == "geneA | geneB"
    assert result.candidate_relation == "5prime_upstream | 3prime_downstream"
    assert result.candidate_distance_bp == "52 | 12"
    assert result.candidate_priority == "high | context_only"
    assert result.candidate_score == "40 | 0"


def test_operon_and_tn5_profile_adjust_score_without_changing_feature_type() -> None:
    gene_a = _gene("geneA", 100, 200, "-", "geneA")
    gene_b = _gene("geneB", 264, 400, "-", "geneB")
    operon_map = {"LT_geneA": "OP001", "LT_geneB": "OP001"}
    result = _annotate(
        252,
        [gene_a, gene_b],
        operon_map=operon_map,
        profile="bidirectional_promoters_and_terminator",
    )

    assert result.feature_relation == "intergenic"
    assert result.candidate_gene_code.startswith("geneA")
    assert result.candidate_score == "80 | 0"
    assert result.candidate_operon_relation == (
        "same_operon_downstream | same_operon_upstream_context"
    )
    assert "terminator_polar_effect" in result.candidate_evidence
    assert result.tn5_structure_profile == "bidirectional_promoters_and_terminator"


def test_500_bp_is_a_reporting_window_not_a_candidate_beyond_it() -> None:
    left = _gene("left", 1, 100, "+")
    right = _gene("right", 1102, 1200, "+")
    result = _annotate(601, [left, right])

    assert result.feature_relation == "intergenic"
    assert result.candidate_gene_code == ""
    assert "No flanking gene falls within the 500-bp" in result.note


def test_gene_exactly_500_bp_away_is_still_reported_as_context() -> None:
    left = _gene("left", 1, 101, "+")
    right = _gene("right", 1102, 1200, "+")
    result = _annotate(601, [left, right])

    assert result.candidate_gene_code == "left"
    assert result.candidate_distance_bp == "500"
    assert result.candidate_priority == "context_only"


def test_load_operon_map_accepts_multiple_gene_identifiers(tmp_path: Path) -> None:
    path = tmp_path / "operons.tsv"
    path.write_text(
        "operon_id\tlocus_tag\tgene_name\nOP1\tLT_A\tgeneA\n",
        encoding="utf-8",
    )

    assert load_operon_map(path) == {"LT_A": "OP1", "geneA": "OP1"}


def test_candidate_functional_annotations_preserve_candidate_alignment() -> None:
    left = _gene("left", 100, 200, "-", "left")
    right = _gene("right", 264, 400, "-", "right")
    genes = {"ctg": [left, right]}
    cdss = {
        "ctg": [
            CDSFeature(
                seqid="ctg",
                start=100,
                end=200,
                strand="-",
                gene_name="left",
                locus_tag="LT_left",
                gene_id="ID_left",
                product="left enzyme",
                ontology_terms=["GO:0001"],
                go_tags=["catalytic activity | GO:0002"],
                attributes={},
            ),
            CDSFeature(
                seqid="ctg",
                start=264,
                end=400,
                strand="-",
                gene_name="right",
                locus_tag="LT_right",
                gene_id="ID_right",
                product="right transporter",
                ontology_terms=["GO:0003"],
                go_tags=["transport | GO:0004"],
                attributes={},
            ),
        ]
    }

    result = annotate_insertion(
        "ctg",
        252,
        genes,
        cdss,
        build_gene_lookup(genes),
        500,
    )

    assert result.candidate_ontology_term == "GO:0001 || GO:0003"
    assert result.candidate_go_tags == (
        "catalytic activity | GO:0002 || transport | GO:0004"
    )
