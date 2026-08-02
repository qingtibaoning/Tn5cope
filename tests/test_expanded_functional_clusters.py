from __future__ import annotations

from tn5cope.functional_clusters import (
    aggregate_gene_records,
    functional_cluster_analysis,
    select_expanded_gene_events,
    tokenize_product,
)


def _base_row(query_id: str, relation: str = "intergenic") -> dict[str, object]:
    return {
        "Query_ID": query_id,
        "Row_number": int(query_id.removeprefix("Q")),
        "Insertion site": 1000,
        "Reference_seqid": "ctg",
        "Hit_status": "unique_hit",
        "Feature_relation": relation,
        "Gene name": "",
        "Gene code": "",
        "Putative function or description": "",
        "Ontology_term": "",
        "Go_tags": "",
        "Candidate gene name": "",
        "Candidate gene code": "",
        "Candidate relation": "",
        "Candidate distance (bp)": "",
        "Candidate strand": "",
        "Candidate priority": "",
        "Candidate function": "",
        "Candidate ontology term": "",
        "Candidate GO tags": "",
        "Candidate score": "",
        "Candidate evidence": "",
        "Candidate operon ID": "",
        "Candidate operon relation": "",
    }


def _candidate_row(
    query_id: str,
    *,
    codes: str,
    priorities: str,
    scores: str,
    products: str = "",
    ontology: str = "",
    go_tags: str = "",
) -> dict[str, object]:
    row = _base_row(query_id)
    count = len(codes.split(" | "))
    row.update(
        {
            "Candidate gene name": "",
            "Candidate gene code": codes,
            "Candidate relation": " | ".join(["5prime_upstream"] * count),
            "Candidate distance (bp)": " | ".join(
                str(25 + index * 100) for index in range(count)
            ),
            "Candidate strand": " | ".join(["+"] * count),
            "Candidate priority": priorities,
            "Candidate function": products or " | ".join(["enzyme"] * count),
            "Candidate ontology term": ontology,
            "Candidate GO tags": go_tags,
            "Candidate score": scores,
            "Candidate evidence": " | ".join(["distance_tier"] * count),
            "Candidate operon ID": " | ".join(["NA"] * count),
            "Candidate operon relation": " | ".join(["not_assessed"] * count),
        }
    )
    return row


def test_unique_top_high_candidate_is_selected() -> None:
    row = _candidate_row(
        "Q001",
        codes="geneA | geneB",
        priorities="high | moderate",
        scores="40 | 25",
        products="transferase | transporter",
        ontology="GO:0001 || GO:0002",
        go_tags="catalytic activity | GO:0003 || transport | GO:0004",
    )

    events, exclusions = select_expanded_gene_events([row])

    assert not exclusions
    assert events[0]["Selected gene ID"] == "geneA"
    assert events[0]["Gene source"] == "candidate_high"
    assert events[0]["Ontology_term"] == "GO:0001"
    assert events[0]["Go_tags"] == "catalytic activity | GO:0003"


def test_unique_top_moderate_candidate_is_selected() -> None:
    row = _candidate_row(
        "Q002",
        codes="geneA | geneB",
        priorities="moderate | context_only",
        scores="25 | 0",
    )

    events, exclusions = select_expanded_gene_events([row])

    assert not exclusions
    assert events[0]["Selected gene ID"] == "geneA"
    assert events[0]["Gene source"] == "candidate_moderate"


def test_tied_top_candidates_are_excluded() -> None:
    row = _candidate_row(
        "Q003",
        codes="geneA | geneB",
        priorities="high | high",
        scores="40 | 40",
    )

    events, exclusions = select_expanded_gene_events([row])

    assert not events
    assert exclusions[0]["Exclusion reason"] == "top_score_tied"


def test_low_priority_top_candidate_is_excluded() -> None:
    row = _candidate_row(
        "Q004",
        codes="geneA",
        priorities="low",
        scores="10",
    )

    events, exclusions = select_expanded_gene_events([row])

    assert not events
    assert exclusions[0]["Exclusion reason"] == "top_priority_not_eligible"


def test_misaligned_candidate_columns_are_excluded() -> None:
    row = _candidate_row(
        "Q005",
        codes="geneA | geneB",
        priorities="high",
        scores="40 | 25",
    )

    events, exclusions = select_expanded_gene_events([row])

    assert not events
    assert exclusions[0]["Exclusion reason"] == "candidate_columns_misaligned"


def test_direct_and_candidate_support_aggregate_to_one_equal_weight_gene() -> None:
    direct = _base_row("Q006", relation="CDS")
    direct.update(
        {
            "Gene name": "abc",
            "Gene code": "geneA",
            "Putative function or description": "ABC transporter ATP-binding protein",
            "Ontology_term": "GO:0005524",
            "Go_tags": "ATP binding | GO:0005524",
        }
    )
    candidate = _candidate_row(
        "Q007",
        codes="geneA",
        priorities="high",
        scores="40",
        products="ABC transporter ATP-binding protein",
        ontology="GO:0005524",
        go_tags="ATP binding | GO:0005524",
    )

    events, exclusions = select_expanded_gene_events([direct, candidate])
    records = aggregate_gene_records(events)
    analysis = functional_cluster_analysis(records)

    assert not exclusions
    assert len(records) == 1
    assignment = analysis["assignments"][0]
    assert assignment["gene_sources"] == "direct; candidate_high"
    assert assignment["source_label"] == "D+H"
    assert assignment["supporting_events"] == 2
    assert assignment["clustering_weight"] == 1


def test_repeated_queries_for_one_gene_remain_one_clustering_unit() -> None:
    rows = []
    for query_id in ("Q008", "Q009"):
        row = _base_row(query_id, relation="CDS")
        row.update(
            {
                "Gene code": "geneA",
                "Putative function or description": "DNA helicase",
            }
        )
        rows.append(row)

    events, _ = select_expanded_gene_events(rows)
    records = aggregate_gene_records(events)

    assert len(events) == 2
    assert len(records) == 1
    assert len(records[0].events) == 2


def test_generic_product_words_do_not_create_functional_similarity() -> None:
    assert tokenize_product("ATP-dependent type IV secretion system protein") == [
        "atp",
        "secretion",
    ]
