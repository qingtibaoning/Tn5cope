from pathlib import Path

from tn5cope.input_parsers import parse_sequences_from_text
from tn5cope.pipeline import load_queries


def test_parse_fasta_text() -> None:
    text = ">a\nacgt\n>b\nnnnn\n"
    assert parse_sequences_from_text(text, "reads.fasta") == ["ACGT", "NNNN"]


def test_parse_csv_text_skips_header() -> None:
    text = "sequence\nacgt\nTGCA\n"
    assert parse_sequences_from_text(text, "reads.csv") == ["ACGT", "TGCA"]


def test_parse_plain_lines() -> None:
    text = "acgt\n\ntgca\n"
    assert parse_sequences_from_text(text, "reads.txt") == ["ACGT", "TGCA"]


def test_cli_query_csv_accepts_header(tmp_path: Path) -> None:
    path = tmp_path / "query.csv"
    path.write_text("sequence\nACGTACGT\n", encoding="ascii")
    records = load_queries(path)
    assert [record.cleaned_sequence for record in records] == ["ACGTACGT"]


def test_cli_query_csv_accepts_no_header(tmp_path: Path) -> None:
    path = tmp_path / "query.csv"
    path.write_text("ACGTACGT\nTGCATGCA\n", encoding="ascii")
    records = load_queries(path)
    assert [record.cleaned_sequence for record in records] == ["ACGTACGT", "TGCATGCA"]
