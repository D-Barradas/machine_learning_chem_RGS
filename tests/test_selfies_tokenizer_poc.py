"""Unit tests for `src/selfies_tokenizer_poc.py`.

These tests are intentionally lightweight and operate on synthetic CSV data.
They verify:
  - The vocab/tokenization utilities behave as expected.
  - The CLI script runs end-to-end on a tiny dataset and produces a checkpoint.

Skip conditions:
  - Tests that require torch will be skipped automatically if torch is not installed.

Run:
  pytest -q

Note: Training is minimal (epochs=1) to keep runtime negligible.
"""

from __future__ import annotations

import csv
import os
import sys
import json
import tempfile
import subprocess
from pathlib import Path

import pytest

# Skip tests if torch not available for CLI run
torch = pytest.importorskip("torch")  # noqa: E305

from src.selfies_tokenizer_poc import (
    split_selfies,
    build_vocab,
    encode_selfies,
)


@pytest.fixture()
def tiny_csvs(tmp_path: Path):
    """Create minimal train/val/test CSVs with SELFIES and DV_C columns."""
    rows = [
        {"text": "[C][O][H]", "DV_C": 0.1},
        {"text": "[C][=O][O][H]", "DV_C": 0.2},
        {"text": "[C][Branch1][=O]", "DV_C": 0.3},
    ]
    splits = {}
    for name in ["train", "val", "test"]:
        file_path = tmp_path / f"hf_dataset_small_{name}.csv"
        with file_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["text", "DV_C"])
            writer.writeheader()
            # Use first 2 for train, 1 for val, 1 for test (reusing rows)
            if name == "train":
                writer.writerow(rows[0])
                writer.writerow(rows[1])
            else:
                writer.writerow(rows[2])
        splits[name] = file_path
    return splits


def test_split_selfies_basic():
    tokens = split_selfies("[C][O][=O]")
    assert tokens == ["[C]", "[O]", "[=O]"]


def test_build_vocab_and_encode():
    import pandas as pd

    df = pd.DataFrame({"text": ["[C][O]", "[C][=O][O]"]})
    stoi, itos = build_vocab(df["text"])
    assert "[C]" in stoi and "[O]" in stoi
    encoded = encode_selfies("[C][O]", stoi, max_len=5)
    # length matches max_len
    assert len(encoded) == 5
    # last two should correspond to tokens (since left-padded)
    assert encoded[-2] == stoi["[C]"] and encoded[-1] == stoi["[O]"]


@pytest.mark.integration
def test_cli_runs_and_creates_checkpoint(tiny_csvs, tmp_path: Path):
    """Run the CLI with epochs=1 and verify checkpoint creation."""
    out_file = tmp_path / "poc_selfies.pth"
    cmd = [
        sys.executable,
        "src/selfies_tokenizer_poc.py",
        "--train",
        str(tiny_csvs["train"]),
        "--val",
        str(tiny_csvs["val"]),
        "--test",
        str(tiny_csvs["test"]),
        "--target",
        "DV_C",
        "--out",
        str(out_file),
        "--epochs",
        "1",
        "--batch_size",
        "2",
        "--embed_dim",
        "16",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"CLI failed: {result.stderr}\n{result.stdout}"
    assert out_file.exists(), "Checkpoint file not created"
    ckpt = torch.load(out_file)
    for key in ["model_state_dict", "stoi", "itos", "args"]:
        assert key in ckpt
    assert isinstance(ckpt["stoi"], dict)
    assert isinstance(ckpt["itos"], list)


def test_checkpoint_vocab_consistency(tiny_csvs, tmp_path: Path):
    """Ensure vocab special tokens exist after a run."""
    out_file = tmp_path / "poc_selfies.pth"
    cmd = [
        sys.executable,
        "src/selfies_tokenizer_poc.py",
        "--train",
        str(tiny_csvs["train"]),
        "--val",
        str(tiny_csvs["val"]),
        "--test",
        str(tiny_csvs["test"]),
        "--target",
        "DV_C",
        "--out",
        str(out_file),
        "--epochs",
        "1",
        "--batch_size",
        "2",
    ]
    subprocess.run(cmd, check=True)
    ckpt = torch.load(out_file)
    stoi = ckpt["stoi"]
    assert "<PAD>" in stoi and "<UNK>" in stoi
