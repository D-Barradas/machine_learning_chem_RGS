#!/usr/bin/env python3
"""
SELFIES tokenizer and vocabulary for RL-based molecule generation.

Builds a vocabulary of SELFIES tokens from a dataset and provides encode/decode.

Usage:
    from selfies_tokenizer import SELFIESTokenizer
    tokenizer = SELFIESTokenizer.from_smiles_file("data/phen_derivatives_scored.csv")
    tokens = tokenizer.encode_selfies(selfies_str)
    selfies_str = tokenizer.decode_tokens(tokens)
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path
from typing import List, Optional

try:
    import selfies as sf
except ImportError:
    sf = None

from rdkit import Chem


class SELFIESTokenizer:
    """Tokenizer for SELFIES strings using bracket notation."""

    def __init__(self, vocab: List[str], pad_token: str = "<PAD>", sos_token: str = "<SOS>", eos_token: str = "<EOS>"):
        self.vocab = vocab
        self.pad_token = pad_token
        self.sos_token = sos_token
        self.eos_token = eos_token

        special = [pad_token, sos_token, eos_token]
        self.all_tokens = special + [t for t in vocab if t not in special]
        self.token_to_id = {t: i for i, t in enumerate(self.all_tokens)}
        self.id_to_token = {i: t for i, t in enumerate(self.all_tokens)}

        self.pad_id = self.token_to_id[pad_token]
        self.sos_id = self.token_to_id[sos_token]
        self.eos_id = self.token_to_id[eos_token]

    @classmethod
    def from_smiles_file(cls, path: str, smiles_col: str = "smiles", max_vocab: int = 200) -> SELFIESTokenizer:
        """Build tokenizer from SMILES in a CSV file."""
        if sf is None:
            raise ImportError("selfies package not installed; run: pip install selfies")

        smiles_list = []
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                s = row.get(smiles_col, "").strip()
                if s:
                    smiles_list.append(s)
                    if len(smiles_list) >= 10000:  # sample for speed
                        break

        # Convert to SELFIES and extract tokens
        token_counter = Counter()
        for smi in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                canonical = Chem.MolToSmiles(mol)
                selfies_str = sf.encoder(canonical)
                tokens = cls._tokenize_selfies(selfies_str)
                token_counter.update(tokens)
            except Exception:
                continue

        # Build vocab from most common tokens
        vocab = [t for t, _ in token_counter.most_common(max_vocab)]
        return cls(vocab)

    @staticmethod
    def _tokenize_selfies(selfies_str: str) -> List[str]:
        """Split SELFIES into bracket tokens."""
        pattern = r"\[[^\]]+\]"
        return re.findall(pattern, selfies_str)

    def encode_selfies(self, selfies_str: str, add_sos: bool = True, add_eos: bool = True) -> List[int]:
        """Encode SELFIES string to token IDs."""
        tokens = self._tokenize_selfies(selfies_str)
        ids = []
        if add_sos:
            ids.append(self.sos_id)
        for t in tokens:
            ids.append(self.token_to_id.get(t, self.token_to_id.get("<PAD>")))
        if add_eos:
            ids.append(self.eos_id)
        return ids

    def decode_tokens(self, token_ids: List[int], remove_special: bool = True) -> str:
        """Decode token IDs to SELFIES string."""
        tokens = [self.id_to_token[i] for i in token_ids if i in self.id_to_token]
        if remove_special:
            tokens = [t for t in tokens if t not in [self.pad_token, self.sos_token, self.eos_token]]
        return "".join(tokens)

    def selfies_to_smiles(self, selfies_str: str) -> Optional[str]:
        """Convert SELFIES to SMILES."""
        if sf is None:
            return None
        try:
            smiles = sf.decoder(selfies_str)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol)
        except Exception:
            return None

    def vocab_size(self) -> int:
        return len(self.all_tokens)
