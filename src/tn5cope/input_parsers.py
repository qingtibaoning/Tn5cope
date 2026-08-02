from __future__ import annotations

import csv
import io
import re
from pathlib import Path
from typing import Iterable, List


def clean_sequence(sequence: str) -> str:
    return re.sub(r"\s+", "", sequence or "").upper()


def parse_fasta_text(text: str) -> List[str]:
    sequences: List[str] = []
    current: List[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                sequences.append(clean_sequence("".join(current)))
                current = []
            continue
        current.append(line)
    if current:
        sequences.append(clean_sequence("".join(current)))
    return [seq for seq in sequences if seq]


def parse_csv_text(text: str) -> List[str]:
    reader = csv.reader(io.StringIO(text))
    rows = [row for row in reader if row]
    if not rows:
        return []
    first_value = clean_sequence(rows[0][0])
    start_index = 1 if first_value and any(base not in "ACGTN" for base in first_value) else 0
    return [clean_sequence(row[0]) for row in rows[start_index:] if clean_sequence(row[0])]


def parse_line_sequences(text: str) -> List[str]:
    return [clean_sequence(line) for line in text.splitlines() if clean_sequence(line)]


def parse_sequences_from_text(text: str, source_name: str = "") -> List[str]:
    stripped = (text or "").strip()
    if not stripped:
        return []
    suffix = Path(source_name).suffix.lower()
    if stripped.startswith(">") or suffix in {".fa", ".fasta", ".fna"}:
        return parse_fasta_text(stripped)
    if suffix == ".csv":
        return parse_csv_text(stripped)
    if "," in stripped.splitlines()[0]:
        return parse_csv_text(stripped)
    return parse_line_sequences(stripped)


def write_query_csv(sequences: Iterable[str], path: Path) -> None:
    seqs = [clean_sequence(sequence) for sequence in sequences if clean_sequence(sequence)]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence"])
        for sequence in seqs:
            writer.writerow([sequence])
