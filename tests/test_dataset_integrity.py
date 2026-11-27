"""Shipped datasets must verify byte-exactly on every platform (SHA-256)."""
import hashlib
import json
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[1]/"skindose_verifier"/"data"
DATASETS = sorted(p.parent for p in DATA.glob("*/metadata.json"))


@pytest.mark.parametrize("directory", DATASETS, ids=[p.name for p in DATASETS])
def test_dataset_files_match_declared_sha256(directory):
    metadata = json.loads((directory/"metadata.json").read_text(encoding="utf-8"))
    assert metadata["files_sha256"]
    for name, expected in metadata["files_sha256"].items():
        assert hashlib.sha256((directory/name).read_bytes()).hexdigest() == expected, (
            f"{directory.name}/{name}: los bytes no coinciden; ¿conversión de finales de línea?")


def test_gitattributes_protects_hashed_datasets():
    attributes = (DATA.parents[1]/".gitattributes").read_text(encoding="utf-8")
    assert "skindose_verifier/data/benmakhlouf_2011/** -text" in attributes
    assert "skindose_verifier/data/derived_*/*.csv -text" in attributes
