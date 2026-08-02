from __future__ import annotations

from tn5cope.position_clusters import (
    MIN_VISIBLE_CLUSTER_SPAN_BP,
    circular_distance,
    density_bin_geometry,
    expanded_display_interval,
    load_expanded_site_records,
    load_site_records,
    position_cluster_analysis,
)


def _main_row(query_id: str, site: int) -> dict[str, object]:
    return {
        "Query_ID": query_id,
        "Insertion site": site,
        "Hit_status": "unique_hit",
        "Feature_relation": "CDS",
        "Gene code": f"gene_{query_id}",
    }


def _expanded_event(query_id: str, site: int) -> dict[str, object]:
    return {
        "Query_ID": query_id,
        "Insertion site": site,
        "Hit_status": "unique_hit",
        "Feature_relation": "CDS",
        "Selected gene code": f"gene_{query_id}",
        "Selected gene ID": f"gene_{query_id}",
    }


def test_expanded_scope_is_clustered_independently_from_all_hits() -> None:
    main_rows = [
        _main_row("Q001", 100),
        _main_row("Q002", 150),
        _main_row("Q003", 40_000),
        _main_row("Q004", 40_050),
    ]
    expanded_rows = [
        _expanded_event("Q001", 100),
        _expanded_event("Q002", 150),
        _expanded_event("Q003", 40_000),
    ]

    all_sites = load_site_records(main_rows)
    expanded_sites = load_expanded_site_records(expanded_rows)
    all_analysis = position_cluster_analysis(all_sites, genome_length=100_000)
    expanded_analysis = position_cluster_analysis(expanded_sites, genome_length=100_000)

    assert len(all_sites) == 4
    assert len(all_analysis["summaries"]) == 2
    assert len(expanded_sites) == 3
    assert len(expanded_analysis["summaries"]) == 1
    assert sum(label < 0 for label in expanded_analysis["labels"]) == 1


def test_expanded_records_deduplicate_coordinates_and_count_queries() -> None:
    rows = [
        _expanded_event("Q001", 500),
        _expanded_event("Q002", 500),
        _expanded_event("Q003", 900),
    ]

    sites = load_expanded_site_records(rows)

    assert [site.site for site in sites] == [500, 900]
    assert sites[0].support_count == 2
    assert sites[0].query_ids == ["Q001", "Q002"]


def test_circular_distance_wraps_across_the_origin() -> None:
    assert circular_distance(10, 990, 1_000) == 20


def test_short_cluster_span_is_widened_only_for_drawing() -> None:
    summary = {
        "interval_start": 10_000,
        "interval_end": 10_005,
        "interval_span_bp": 5,
    }

    display_start, display_end = expanded_display_interval(summary, 100_000)

    assert display_end - display_start == MIN_VISIBLE_CLUSTER_SPAN_BP
    assert summary["interval_start"] == 10_000
    assert summary["interval_end"] == 10_005


def test_small_genome_density_bin_stays_inside_sector() -> None:
    centers, widths = density_bin_geometry(600)

    assert centers.tolist() == [300.0]
    assert widths.tolist() == [552.0]


def test_small_genome_cluster_display_is_capped_to_genome() -> None:
    summary = {
        "interval_start": 151,
        "interval_end": 459,
        "interval_span_bp": 308,
    }

    assert expanded_display_interval(summary, 600) == (0, 600)
