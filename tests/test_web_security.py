from __future__ import annotations

import pytest


pytest.importorskip("fastapi")

from tn5cope.web import safe_path_component


def test_download_path_components_accept_generated_names() -> None:
    assert safe_path_component("20260803_tail_pcr_a1b2c3d4")
    assert safe_path_component("tn5_results.xlsx")


@pytest.mark.parametrize(
    "value",
    ["", ".", "..", "../README.md", ".hidden", "nested/file.tsv"],
)
def test_download_path_components_reject_traversal(value: str) -> None:
    assert not safe_path_component(value)
