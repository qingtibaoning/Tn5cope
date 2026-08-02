from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, lines as mlines, patches as mpatches
from matplotlib.text import Text
import numpy as np
from pycirclize import Circos
from sklearn.cluster import DBSCAN


VALID_HIT_STATUSES = {"unique_hit", "rescued_unique_hit"}
POSITION_DBSCAN_EPS = 20_000
POSITION_DBSCAN_MIN_SAMPLES = 2
POSITION_SENSITIVITY_EPS = (5_000, 10_000, 20_000)
POSITION_SENSITIVITY_MIN_SAMPLES = (2, 3)
POSITION_DENSITY_BIN_SIZE = 50_000
PLOT_DPI = 220
GENOME_TICK_INTERVAL = 1_000_000
FIGURE_SIZE = (10.8, 11.4)
FONT_FAMILY = "Times New Roman"
MIN_VISIBLE_CLUSTER_SPAN_BP = 7_000
TIMES_NEW_ROMAN_FONT_PATHS = (
    Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
    Path("/Library/Fonts/Times New Roman.ttf"),
    Path.home() / "Library/Fonts/Times New Roman.ttf",
)
TIMES_NEW_ROMAN_BOLD_FONT_PATHS = (
    Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"),
    Path("/Library/Fonts/Times New Roman Bold.ttf"),
    Path.home() / "Library/Fonts/Times New Roman Bold.ttf",
)

PLOT_COLORS = {
    "outer_ring": "#000000",
    "inner_ring": "#b8b8b8",
    "text": "#1c1c1c",
    "muted_text": "#777777",
    "expanded_cluster": "#b8a5f5",
    "expanded_clustered_site": "#4ec57a",
    "expanded_unclustered_site": "#7e8790",
    "density": "#ee8ebb",
    "plus_strand_gene": "#74b8f2",
    "minus_strand_gene": "#f2a36b",
}

ASSIGNMENT_HEADERS = [
    "site",
    "gene_label",
    "feature_relation",
    "support_count",
    "query_ids",
    "cluster_id",
]
SUMMARY_HEADERS = [
    "cluster_id",
    "site_count",
    "total_support",
    "interval_start",
    "interval_end",
    "interval_span_bp",
    "member_sites",
    "member_genes",
]
SENSITIVITY_HEADERS = ["eps_bp", "min_samples", "cluster_count", "noise_site_count"]


@dataclass
class SiteRecord:
    site: int
    feature_relation: str
    gene_label: str
    support_count: int
    query_ids: List[str]


def output_path(output_prefix: Path, suffix: str) -> Path:
    return output_prefix.with_name(output_prefix.name + suffix)


def load_site_records(rows: Sequence[Dict[str, object]]) -> List[SiteRecord]:
    grouped: Dict[int, SiteRecord] = {}
    for row in rows:
        if str(row.get("Hit_status", "")) not in VALID_HIT_STATUSES:
            continue
        raw_site = row.get("Insertion site", "")
        if raw_site in ("", None):
            continue
        site = int(raw_site)
        gene_label = (
            str(row.get("Gene code") or row.get("Gene name") or "").strip()
            or "intergenic"
        )
        grouped.setdefault(
            site,
            SiteRecord(
                site=site,
                feature_relation=str(row.get("Feature_relation") or "unknown"),
                gene_label=gene_label,
                support_count=0,
                query_ids=[],
            ),
        )
        record = grouped[site]
        record.support_count += 1
        query_id = str(row.get("Query_ID") or "")
        if query_id:
            record.query_ids.append(query_id)
    return sorted(grouped.values(), key=lambda item: item.site)


def load_expanded_site_records(
    expanded_event_rows: Sequence[Dict[str, object]],
) -> List[SiteRecord]:
    """Collapse Expanded candidate events to unique insertion coordinates.

    Expanded events have already passed the formal inclusion rule: reliable direct
    gene hits plus a unique highest-scoring high/moderate intergenic candidate.
    Position clustering uses the selected events, not every high-confidence genomic
    alignment, so the PC labels and site ring match the formal reporting scope.
    """

    grouped: Dict[int, SiteRecord] = {}
    for row in expanded_event_rows:
        if str(row.get("Hit_status", "")) not in VALID_HIT_STATUSES:
            continue
        raw_site = row.get("Insertion site", "")
        if raw_site in ("", None):
            continue
        site = int(raw_site)
        gene_label = (
            str(
                row.get("Selected gene code")
                or row.get("Selected gene name")
                or row.get("Selected gene ID")
                or ""
            ).strip()
            or "unresolved"
        )
        grouped.setdefault(
            site,
            SiteRecord(
                site=site,
                feature_relation=str(row.get("Feature_relation") or "unknown"),
                gene_label=gene_label,
                support_count=0,
                query_ids=[],
            ),
        )
        record = grouped[site]
        record.support_count += 1
        query_id = str(row.get("Query_ID") or "")
        if query_id:
            record.query_ids.append(query_id)
    return sorted(grouped.values(), key=lambda item: item.site)


def circular_distance(a: int, b: int, genome_length: int) -> int:
    delta = abs(a - b)
    return min(delta, genome_length - delta)


def build_distance_matrix(sites: Sequence[SiteRecord], genome_length: int) -> np.ndarray:
    n_sites = len(sites)
    matrix = np.zeros((n_sites, n_sites), dtype=float)
    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            distance = circular_distance(sites[i].site, sites[j].site, genome_length)
            matrix[i, j] = distance
            matrix[j, i] = distance
    return matrix


def circular_interval(coords: Sequence[int], genome_length: int) -> Tuple[int, int, int]:
    if not coords:
        return 0, 0, 0
    ordered = sorted(coords)
    if len(ordered) == 1:
        coord = ordered[0]
        return coord, coord, 0
    gaps = []
    for left, right in zip(ordered, ordered[1:]):
        gaps.append((right - left, left, right))
    gaps.append((genome_length - ordered[-1] + ordered[0], ordered[-1], ordered[0]))
    largest_gap, gap_start, gap_end = max(gaps, key=lambda item: item[0])
    return gap_end, gap_start, genome_length - largest_gap


def cluster_labels(distance_matrix: np.ndarray, eps: int, min_samples: int) -> np.ndarray:
    if distance_matrix.shape[0] == 0:
        return np.array([], dtype=int)
    model = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    return model.fit_predict(distance_matrix)


def position_cluster_analysis(sites: Sequence[SiteRecord], genome_length: int) -> Dict[str, object]:
    distance_matrix = build_distance_matrix(sites, genome_length)
    labels = cluster_labels(distance_matrix, POSITION_DBSCAN_EPS, POSITION_DBSCAN_MIN_SAMPLES)

    assignments: List[Dict[str, object]] = []
    cluster_members: Dict[int, List[SiteRecord]] = defaultdict(list)
    for site, label in zip(sites, labels):
        cluster_members[int(label)].append(site)
        assignments.append(
            {
                "site": site.site,
                "gene_label": site.gene_label,
                "feature_relation": site.feature_relation,
                "support_count": site.support_count,
                "query_ids": "; ".join(site.query_ids),
                "cluster_id": f"PC{int(label) + 1:02d}" if label >= 0 else "NOISE",
            }
        )

    summaries: List[Dict[str, object]] = []
    for label, members in sorted(cluster_members.items()):
        if label < 0:
            continue
        coords = [member.site for member in members]
        start, end, span = circular_interval(coords, genome_length)
        summaries.append(
            {
                "cluster_id": f"PC{label + 1:02d}",
                "site_count": len(members),
                "total_support": sum(member.support_count for member in members),
                "interval_start": start,
                "interval_end": end,
                "interval_span_bp": span,
                "member_sites": "; ".join(
                    str(member.site) for member in sorted(members, key=lambda item: item.site)
                ),
                "member_genes": "; ".join(sorted(member.gene_label for member in members)),
            }
        )

    sensitivity_rows: List[Dict[str, object]] = []
    for eps in POSITION_SENSITIVITY_EPS:
        for min_samples in POSITION_SENSITIVITY_MIN_SAMPLES:
            trial_labels = cluster_labels(distance_matrix, eps, min_samples)
            sensitivity_rows.append(
                {
                    "eps_bp": eps,
                    "min_samples": min_samples,
                    "cluster_count": len({label for label in trial_labels if label >= 0}),
                    "noise_site_count": sum(1 for label in trial_labels if label < 0),
                }
            )

    return {
        "assignments": assignments,
        "summaries": summaries,
        "sensitivity": sensitivity_rows,
        "labels": labels,
    }


def write_tsv(path: Path, headers: Sequence[str], rows: Sequence[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(headers), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def setup_times_new_roman() -> str:
    for path in TIMES_NEW_ROMAN_FONT_PATHS + TIMES_NEW_ROMAN_BOLD_FONT_PATHS:
        if path.exists():
            font_manager.fontManager.addfont(str(path))
    for path in TIMES_NEW_ROMAN_FONT_PATHS:
        if path.exists():
            font_name = font_manager.FontProperties(fname=str(path)).get_name()
            break
    else:
        font_name = FONT_FAMILY
    matplotlib.rcParams.update(
        {
            "font.family": font_name,
            "font.serif": [font_name],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )
    return font_name


def cluster_midpoint(summary: Dict[str, object], genome_length: int) -> float:
    start = int(summary["interval_start"])
    end = int(summary["interval_end"])
    if start <= end:
        return (start + end) / 2
    return ((start + end + genome_length) / 2) % genome_length


def cluster_label_radius(summary: Dict[str, object], index: int) -> float:
    manual_radii = {
        "PC01": 103,
        "PC02": 109,
        "PC03": 103,
        "PC04": 109,
        "PC05": 104,
        "PC06": 110,
        "PC07": 104,
        "PC08": 110,
        "PC09": 104,
        "PC10": 110,
        "PC11": 104,
        "PC12": 110,
        "PC13": 104,
        "PC14": 110,
    }
    cluster_id = str(summary["cluster_id"])
    return manual_radii.get(cluster_id, (103, 107, 111)[index % 3])


def add_genome_axis(sector, genome_seqid: str, genome_length: int, font_name: str) -> None:
    track = sector.add_track((98, 100))
    track.axis(fc="none", ec=PLOT_COLORS["outer_ring"], lw=0.9)
    track.xticks_by_interval(
        GENOME_TICK_INTERVAL,
        tick_length=1.7,
        label_size=8.5,
        label_formatter=lambda value: f"{value / 1_000_000:.1f} Mb",
        line_kws={"color": PLOT_COLORS["outer_ring"], "lw": 0.65},
        text_kws={"fontname": font_name, "color": PLOT_COLORS["muted_text"]},
    )


def expanded_display_interval(
    summary: Dict[str, object],
    genome_length: int,
) -> Tuple[int, int]:
    """Widen tiny PC spans only for drawing; tables retain the exact interval."""

    start = int(summary["interval_start"])
    end = int(summary["interval_end"])
    span = int(summary["interval_span_bp"])
    visible_span = min(MIN_VISIBLE_CLUSTER_SPAN_BP, genome_length)
    if span >= visible_span:
        return start, end
    if visible_span == genome_length:
        return 0, genome_length
    midpoint = cluster_midpoint(summary, genome_length)
    display_start = int(round(midpoint - visible_span / 2))
    display_end = int(round(midpoint + visible_span / 2))
    return display_start % genome_length, display_end % genome_length


def draw_cluster_span(
    track,
    summary: Dict[str, object],
    r_lim: Tuple[float, float],
    genome_length: int,
) -> None:
    start, end = expanded_display_interval(summary, genome_length)
    if start <= end:
        track.rect(
            start,
            end,
            r_lim=r_lim,
            fc=PLOT_COLORS["expanded_cluster"],
            ec="none",
            alpha=0.58,
        )
    else:
        track.rect(
            start,
            genome_length,
            r_lim=r_lim,
            fc=PLOT_COLORS["expanded_cluster"],
            ec="none",
            alpha=0.58,
        )
        track.rect(
            0,
            end,
            r_lim=r_lim,
            fc=PLOT_COLORS["expanded_cluster"],
            ec="none",
            alpha=0.58,
        )


def density_bin_geometry(genome_length: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return in-range bar centers and widths, including a shortened final bin."""

    if genome_length <= 0:
        raise ValueError("genome_length must be positive")
    starts = np.arange(0, genome_length, POSITION_DENSITY_BIN_SIZE, dtype=float)
    ends = np.minimum(starts + POSITION_DENSITY_BIN_SIZE, genome_length)
    centers = (starts + ends) / 2
    widths = (ends - starts) * 0.92
    return centers, widths


def add_density_track(
    sector,
    sites: Sequence[SiteRecord],
    genome_length: int,
    font_name: str,
) -> int:
    track = sector.add_track((26, 46))
    track.axis(fc="none", ec=PLOT_COLORS["inner_ring"], lw=0.5)
    x, widths = density_bin_geometry(genome_length)
    n_bins = len(x)
    counts = np.zeros(n_bins)
    for site in sites:
        counts[(site.site - 1) // POSITION_DENSITY_BIN_SIZE] += 1
    vmax = float(max(counts.max(), 1))
    track.bar(
        x,
        counts,
        width=widths,
        vmin=0,
        vmax=vmax,
        color=PLOT_COLORS["density"],
        ec="none",
        alpha=0.68,
    )
    track.yticks(
        [0, vmax],
        labels=["0", str(int(vmax))],
        vmin=0,
        vmax=vmax,
        label_size=7.5,
        text_kws={"fontname": font_name, "color": PLOT_COLORS["muted_text"]},
    )
    return int(vmax)


def add_gene_tracks(sector, genes: Sequence[object]) -> None:
    plus_track = sector.add_track((55, 59))
    minus_track = sector.add_track((49, 53))
    for track in (plus_track, minus_track):
        track.axis(fc="none", ec=PLOT_COLORS["inner_ring"], lw=0.35)

    for gene in genes:
        start = int(getattr(gene, "start"))
        end = int(getattr(gene, "end"))
        strand = str(getattr(gene, "strand", ""))
        if strand == "+":
            plus_track.rect(
                start,
                end,
                r_lim=(55.4, 58.6),
                fc=PLOT_COLORS["plus_strand_gene"],
                ec="none",
                alpha=0.64,
            )
        elif strand == "-":
            minus_track.rect(
                start,
                end,
                r_lim=(49.4, 52.6),
                fc=PLOT_COLORS["minus_strand_gene"],
                ec="none",
                alpha=0.64,
            )


def add_circle_legend(
    fig,
    font_name: str,
    cluster_count: int,
    density_max: int,
) -> None:
    cluster_range = f"PC01-PC{cluster_count:02d}" if cluster_count else "none"
    handles = [
        mlines.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markerfacecolor=PLOT_COLORS["expanded_clustered_site"],
            markeredgecolor="white",
            color="none",
            label="Expanded candidate - clustered",
            markersize=7,
        ),
        mlines.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markerfacecolor="white",
            markeredgecolor=PLOT_COLORS["expanded_unclustered_site"],
            color="none",
            label="Expanded candidate - unclustered",
            markersize=7,
        ),
        mpatches.Patch(
            color=PLOT_COLORS["expanded_cluster"],
            alpha=0.58,
            label=f"Expanded candidate clusters ({cluster_range})",
        ),
        mpatches.Patch(
            color=PLOT_COLORS["density"],
            alpha=0.68,
            label=(
                "50 kb density: all unique high-confidence loci "
                f"(maximum {density_max})"
            ),
        ),
        mlines.Line2D(
            [],
            [],
            color=PLOT_COLORS["plus_strand_gene"],
            lw=4,
            label="+ strand genes",
        ),
        mlines.Line2D(
            [],
            [],
            color=PLOT_COLORS["minus_strand_gene"],
            lw=4,
            label="- strand genes",
        ),
    ]
    legend = fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.060),
        ncol=2,
        frameon=False,
        prop={"family": font_name, "size": 9.5},
        handlelength=1.7,
        columnspacing=1.4,
        labelspacing=0.8,
    )
    for text in legend.get_texts():
        text.set_fontname(font_name)


def apply_font(fig, font_name: str) -> None:
    for text in fig.findobj(match=Text):
        text.set_fontname(font_name)


def plot_circle(
    all_high_confidence_sites: Sequence[SiteRecord],
    expanded_sites: Sequence[SiteRecord],
    expanded_position_analysis: Dict[str, object],
    genome_length: int,
    output_file: Path,
    genome_seqid: str = "Reference genome",
    genes: Sequence[object] = (),
) -> None:
    font_name = setup_times_new_roman()
    circos = Circos({genome_seqid: genome_length}, start=-270, end=90, space=0)
    sector = circos.sectors[0]
    add_genome_axis(sector, genome_seqid, genome_length, font_name)

    if not expanded_sites:
        fig = circos.plotfig(dpi=PLOT_DPI, figsize=FIGURE_SIZE)
        fig.suptitle(
            "Expanded-candidate Tn5 insertion circle",
            y=0.988,
            fontsize=18,
            fontname=font_name,
            fontweight="bold",
        )
        fig.text(
            0.5,
            0.5,
            "No Expanded candidate insertion loci.",
            ha="center",
            va="center",
            fontname=font_name,
        )
        add_circle_legend(fig, font_name, 0, 0)
        apply_font(fig, font_name)
        fig.savefig(output_file, dpi=300, bbox_inches="tight", pad_inches=0.35)
        fig.savefig(output_file.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.35)
        fig.savefig(output_file.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.35)
        plt.close(fig)
        return

    summaries = list(expanded_position_analysis["summaries"])
    cluster_track = sector.add_track((76, 94))
    cluster_track.axis(fc="none", ec=PLOT_COLORS["inner_ring"], lw=0.5)
    for summary in summaries:
        draw_cluster_span(
            cluster_track,
            summary,
            (79, 92),
            genome_length,
        )
    for idx, summary in enumerate(summaries):
        cluster_track.text(
            f"{summary['cluster_id']} ({summary['site_count']} hits)",
            x=cluster_midpoint(summary, genome_length),
            r=cluster_label_radius(summary, idx),
            size=8.5,
            color=PLOT_COLORS["text"],
            fontname=font_name,
        )

    point_track = sector.add_track((63, 72))
    point_track.axis(fc="none", ec=PLOT_COLORS["inner_ring"], lw=0.5)
    labels = list(expanded_position_analysis["labels"])
    clustered_sites = [
        site for site, label in zip(expanded_sites, labels) if int(label) >= 0
    ]
    unclustered_sites = [
        site for site, label in zip(expanded_sites, labels) if int(label) < 0
    ]
    if clustered_sites:
        point_track.scatter(
            [site.site for site in clustered_sites],
            [0.5] * len(clustered_sites),
            vmin=0,
            vmax=1,
            s=[22 + 16 * site.support_count for site in clustered_sites],
            color=PLOT_COLORS["expanded_clustered_site"],
            ec="white",
            lw=0.45,
            alpha=0.95,
        )
    if unclustered_sites:
        point_track.scatter(
            [site.site for site in unclustered_sites],
            [0.5] * len(unclustered_sites),
            vmin=0,
            vmax=1,
            s=[22 + 16 * site.support_count for site in unclustered_sites],
            color="white",
            ec=PLOT_COLORS["expanded_unclustered_site"],
            lw=0.75,
            alpha=0.95,
        )

    add_gene_tracks(sector, genes)
    density_max = add_density_track(
        sector,
        all_high_confidence_sites,
        genome_length,
        font_name,
    )

    fig = circos.plotfig(dpi=PLOT_DPI, figsize=FIGURE_SIZE)
    fig.suptitle(
        "Expanded-candidate Tn5 insertion circle",
        y=0.988,
        fontsize=18,
        fontname=font_name,
        fontweight="bold",
    )
    query_count = sum(site.support_count for site in expanded_sites)
    noise_count = sum(1 for label in labels if int(label) < 0)
    fig.text(
        0.5,
        0.535,
        genome_seqid,
        ha="center",
        va="center",
        fontname=font_name,
        fontsize=11,
        fontweight="bold",
        color=PLOT_COLORS["text"],
    )
    fig.text(
        0.5,
        0.511,
        "Expanded candidate scope",
        ha="center",
        va="center",
        fontname=font_name,
        fontsize=10,
        color=PLOT_COLORS["text"],
    )
    fig.text(
        0.5,
        0.487,
        f"{len(expanded_sites)} unique sites - {query_count} queries",
        ha="center",
        va="center",
        fontname=font_name,
        fontsize=10,
        color=PLOT_COLORS["text"],
    )
    fig.text(
        0.5,
        0.463,
        f"{len(summaries)} position clusters - {noise_count} unclustered",
        ha="center",
        va="center",
        fontname=font_name,
        fontsize=9.5,
        color=PLOT_COLORS["muted_text"],
    )
    add_circle_legend(fig, font_name, len(summaries), density_max)
    apply_font(fig, font_name)
    fig.savefig(output_file, dpi=300, bbox_inches="tight", pad_inches=0.35)
    fig.savefig(output_file.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.35)
    fig.savefig(output_file.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)


def write_methods_file(
    output_prefix: Path,
    genome_seqid: str,
    genome_length: int,
    all_site_count: int,
    expanded_site_count: int,
    expanded_query_count: int,
) -> None:
    text = f"""Genome position clustering

Two reported scopes
1. All high-confidence position clustering
- Source: {output_path(output_prefix, "_main_results.tsv").name}
- Included hit_status: unique_hit, rescued_unique_hit
- Deduplication unit: unique insertion coordinate
- Unique loci: {all_site_count}
- Outputs: *_position_cluster_assignments.tsv, *_position_cluster_summary.tsv,
  *_position_cluster_sensitivity.tsv

2. Expanded candidate position clustering
- Source: {output_path(output_prefix, "_expanded_gene_events.tsv").name}
- Included events: reliable direct gene hits plus uniquely top-ranked high/moderate
  intergenic candidates
- Deduplication unit: unique insertion coordinate
- Unique loci: {expanded_site_count}
- Supporting query events: {expanded_query_count}
- Outputs: *_expanded_position_cluster_assignments.tsv,
  *_expanded_position_cluster_summary.tsv,
  *_expanded_position_cluster_sensitivity.tsv

Reference genome:
- seqid = {genome_seqid}
- genome_length_bp = {genome_length}

Distance definition:
- circular genomic distance d(a,b) = min(|a-b|, genome_length - |a-b|)

Primary clustering:
- algorithm = DBSCAN
- metric = precomputed circular distance matrix
- eps_bp = {POSITION_DBSCAN_EPS}
- min_samples = {POSITION_DBSCAN_MIN_SAMPLES}

Sensitivity grid:
- eps_bp tested = {", ".join(str(value) for value in POSITION_SENSITIVITY_EPS)}
- min_samples tested = {", ".join(str(value) for value in POSITION_SENSITIVITY_MIN_SAMPLES)}

Circle plot layers:
- renderer = pyCirclize Circos-style circular plot
- all text = Times New Roman when available
- outer double boundary = black
- inner structural rings = gray
- outer PC spans and PC labels = Expanded candidate position clusters only
- PC label format = PCxx (n hits), no leader lines
- point ring = Expanded candidate loci; filled green for clustered, hollow gray for
  unclustered; point size scales with support_count
- gene tracks = plus strand blue, minus strand orange
- inner pink bar ring = all unique high-confidence loci per
  {POSITION_DENSITY_BIN_SIZE} bp bin
- short PC intervals are widened to {MIN_VISIBLE_CLUSTER_SPAN_BP} bp only for
  visibility; exact coordinates remain in the TSV files
- image outputs = PNG, PDF, and SVG
"""
    output_path(output_prefix, "_position_clustering_methods.txt").write_text(text, encoding="utf-8")


def run_position_cluster_analysis(
    main_rows: Sequence[Dict[str, object]],
    expanded_event_rows: Sequence[Dict[str, object]],
    genome_seqid: str,
    genome_length: int,
    output_prefix: Path,
    genes: Sequence[object] = (),
) -> Dict[str, object]:
    all_sites = load_site_records(main_rows)
    all_analysis = position_cluster_analysis(all_sites, genome_length)
    expanded_sites = load_expanded_site_records(expanded_event_rows)
    expanded_analysis = position_cluster_analysis(expanded_sites, genome_length)

    write_tsv(
        output_path(output_prefix, "_position_cluster_assignments.tsv"),
        ASSIGNMENT_HEADERS,
        all_analysis["assignments"],
    )
    write_tsv(
        output_path(output_prefix, "_position_cluster_summary.tsv"),
        SUMMARY_HEADERS,
        all_analysis["summaries"],
    )
    write_tsv(
        output_path(output_prefix, "_position_cluster_sensitivity.tsv"),
        SENSITIVITY_HEADERS,
        all_analysis["sensitivity"],
    )
    write_tsv(
        output_path(output_prefix, "_expanded_position_cluster_assignments.tsv"),
        ASSIGNMENT_HEADERS,
        expanded_analysis["assignments"],
    )
    write_tsv(
        output_path(output_prefix, "_expanded_position_cluster_summary.tsv"),
        SUMMARY_HEADERS,
        expanded_analysis["summaries"],
    )
    write_tsv(
        output_path(output_prefix, "_expanded_position_cluster_sensitivity.tsv"),
        SENSITIVITY_HEADERS,
        expanded_analysis["sensitivity"],
    )
    plot_circle(
        all_sites,
        expanded_sites,
        expanded_analysis,
        genome_length,
        output_path(output_prefix, "_genome_circle_plot.png"),
        genome_seqid,
        genes,
    )
    write_methods_file(
        output_prefix,
        genome_seqid,
        genome_length,
        len(all_sites),
        len(expanded_sites),
        sum(site.support_count for site in expanded_sites),
    )

    all_labels = all_analysis["labels"]
    expanded_labels = expanded_analysis["labels"]
    return {
        "unique_insertion_loci": len(all_sites),
        "position_cluster_count": len(all_analysis["summaries"]),
        "position_noise_loci": sum(1 for label in all_labels if label < 0),
        "expanded_unique_insertion_loci": len(expanded_sites),
        "expanded_position_cluster_count": len(expanded_analysis["summaries"]),
        "expanded_position_noise_loci": sum(
            1 for label in expanded_labels if label < 0
        ),
        "all_assignment_rows": list(all_analysis["assignments"]),
        "all_summary_rows": list(all_analysis["summaries"]),
        "expanded_assignment_rows": list(expanded_analysis["assignments"]),
        "expanded_summary_rows": list(expanded_analysis["summaries"]),
    }
