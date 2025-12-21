#!/usr/bin/env python3
"""
Generate orthophenanthroline derivatives (positions 2–9) with specified substituents
and emit SMILES plus metadata to a CSV.

- Core: 1,10-phenanthroline (PubChem CID 1318)
- Substitution sites: aromatic carbons bearing H (8 sites ≈ positions 2–9)
- Substituents: user-provided names mapped to attachment SMILES with a dummy '[*]'
- Output: CSV with columns: system, smiles, sites, substituents, base_smiles, n_subs, dv_c_sum (optional), notes

Notes
- Atom numbering 2–9 is approximated by ordering detected aromatic C-H sites; this is stable
  but may differ from IUPAC numbering. We include RDKit atom indices for reproducibility.
- Attachment SMILES use a leading '[*]' to indicate the atom that bonds to the ring carbon.
- Requires RDKit.

Example
    python src/generate_phen_derivatives.py \
        --out data/phen_derivatives_generated.csv \
        --subs NO2 CN CHO OCF3 Cl Br NH2 OMe Me Et SiH3 \
        --modes mono di

"""
from __future__ import annotations
import argparse
import csv
import itertools
import os
import random
from typing import Dict, List, Tuple, Optional

from rdkit import Chem
from rdkit.Chem import AllChem

BASE_SMILES = "C1=CC2=C(C3=C(C=CC=N3)C=C2)N=C1"  # 1,10-phenanthroline

# Map common substituent names to attachment SMILES (with '[*]' dummy as attachment point)
SUB_ATTACH_SMILES: Dict[str, str] = {
    # Electron-withdrawing
    "NO2": "[*][N+](=O)[O-]",  # nitro via N (charged nitro)
    "CN": "[*]C#N",          # cyano via carbon
    "CHO": "[*]C=O",         # formyl via carbonyl carbon
    "OCF3": "[*]OC(F)(F)F",  # O-trifluoromethyl via O
    "Cl": "[*]Cl",
    "Br": "[*]Br",
    # Electron-donating
    "NH2": "[*]N",           # amino via N (implicit Hs)
    "OMe": "[*]OC",          # methoxy via O
    "Me": "[*]C",            # methyl via carbon
    "Et": "[*]CC",           # ethyl via carbon
    "SiH3": "[*][SiH3]",     # silicon with explicit hydrogens
}

# Optional DV_C values for a subset (from Table3 in e_int_consolidated.csv if available)
DEFAULT_DV_C: Dict[str, float] = {
    # Only for known entries; extendable by parsing e_int_consolidated.csv
    "CHO": 13.8,
    "CN": 18.0,
    "OCN": 14.2,
    "OCOCF3": 12.2,
    "ONO2": 13.4,
    "CF2CF3": 14.7,
    "Cyclo-C4F7": 16.1,
    "SO2Me": 16.8,
    "OMe": 5.0,
    "Et": 4.0,
    "Cyclohexyl": 3.8,
    "NH2": 9.0,
    "NHNH2": 11.2,
    "NMe2": 12.3,
    "O(CH2)3CH3": 5.9,
    "OCHMe2": 6.5,
    "SiClMe2": 4.8,
    "SCHQCH2": 4.1,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate orthophenanthroline derivatives")
    p.add_argument("--out", required=True, help="Output CSV path")
    p.add_argument("--subs", nargs="+", help="List of substituent names", default=list(SUB_ATTACH_SMILES.keys()))
    p.add_argument("--modes", nargs="+", choices=["mono", "di", "tri", "tetra"], default=["mono", "di"],
                   help="Substitution multiplicities to generate")
    p.add_argument("--eint", help="Optional path to e_int_consolidated.csv to enrich DV_C mapping")
    p.add_argument(
        "--identical-only",
        action="store_true",
        help="Generate only identical-substituent combinations for a given multiplicity to limit explosion",
    )
    p.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional cap on total generated rows; stops early once reached",
    )
    p.add_argument(
        "--shuffle-seed",
        type=int,
        default=13,
        help="Shuffle seed to randomize position/substituent order when capping",
    )
    return p.parse_args()


def get_substituent_smiles(sub_name: str) -> str:
    if sub_name not in SUB_ATTACH_SMILES:
        raise ValueError(f"No attachment SMILES defined for substituent '{sub_name}'")
    return SUB_ATTACH_SMILES[sub_name]


def load_dv_c_from_eint(path: str) -> Dict[str, float]:
    dv_map: Dict[str, float] = {}
    if not path or not os.path.exists(path):
        return dv_map
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("Source", "") == "Table3_eq8":
                name = row.get("System", "").strip()
                try:
                    dv = float(row.get("DV_C (kcal/mol)", "").strip())
                except Exception:
                    dv = None
                if dv is not None and name:
                    dv_map[name] = dv
            # Also allow bare substituent names that match keys in DEFAULT_DV_C
            # The 'System' column is a substituent label for Table3 rows.
    return dv_map


def load_dv_c_from_mapping(path: str = "data/substituent_dv_mappings.csv") -> Dict[str, float]:
    dv_map: Dict[str, float] = {}
    try:
        if not os.path.exists(path):
            return dv_map
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = (row.get("substituent") or "").strip()
                dv_s = (row.get("DV_C") or "").strip()
                if not name or not dv_s:
                    continue
                try:
                    dv = float(dv_s)
                except Exception:
                    continue
                dv_map[name] = dv
    except Exception:
        # If anything goes wrong, silently fall back to defaults
        return {}
    return dv_map


def find_aromatic_CH_sites(mol: Chem.Mol) -> List[int]:
    """Return atom indices of aromatic carbons bearing implicit H (C-H sites)."""
    sites = []
    for a in mol.GetAtoms():
        if a.GetIsAromatic() and a.GetAtomicNum() == 6:
            # Approximate: available for substitution if it has one implicit H
            hcount = a.GetTotalNumHs()
            if hcount > 0:
                sites.append(a.GetIdx())
    # Stable order
    return sorted(sites)


def attach_substituent(mol: Chem.Mol, ring_idx: int, attach_smiles: str) -> Chem.Mol:
    """Attach substituent at given aromatic carbon index using '[*]' dummy in attach_smiles."""
    subs = Chem.MolFromSmiles(attach_smiles)
    if subs is None:
        raise ValueError(f"Invalid substituent SMILES: {attach_smiles}")
    # Find the dummy atom
    dummy_idx = None
    for a in subs.GetAtoms():
        if a.GetAtomicNum() == 0:  # '*'
            dummy_idx = a.GetIdx()
            break
    if dummy_idx is None:
        raise ValueError("Attachment SMILES must contain a '[*]' dummy atom")

    combo = Chem.CombineMols(mol, subs)
    em = Chem.EditableMol(combo)
    offset = mol.GetNumAtoms()
    em.AddBond(ring_idx, offset + dummy_idx, Chem.BondType.SINGLE)
    # Remove the dummy atom; need to track index: removing requires updated index
    # To remove safely, we first convert to RWMol and then remove by atom mapping.
    tmp = em.GetMol()
    # Recompute dummy index location (it may have shifted but here it's the same)
    em2 = Chem.EditableMol(tmp)
    em2.RemoveAtom(offset + dummy_idx)
    out = em2.GetMol()
    Chem.SanitizeMol(out)
    return out


def generate_derivatives(
    sub_names: List[str],
    modes: List[str],
    eint_path: str | None = None,
    identical_only: bool = False,
    max_rows: Optional[int] = None,
    shuffle_seed: int = 13,
) -> List[Dict[str, str]]:
    base = Chem.MolFromSmiles(BASE_SMILES)
    if base is None:
        raise RuntimeError("Failed to parse base phenanthroline SMILES")
    sites = find_aromatic_CH_sites(base)
    # Map sites to pseudo positions 2..9
    pos_labels = {idx: str(2 + i) for i, idx in enumerate(sites)}

    # DV_C mapping
    dv_from_eint = load_dv_c_from_eint(eint_path) if eint_path else {}
    dv_from_mapping = load_dv_c_from_mapping()

    rows: List[Dict[str, str]] = []
    # Preload substituent smiles
    sub_attach = {name: get_substituent_smiles(name) for name in sub_names}

    rng = random.Random(shuffle_seed)

    # Helper to label system string like '5-NO2-phen' or '3,8-(OMe)2-phen'
    def system_label(pos_list: List[int], subs_list: List[str]) -> str:
        # If all substituents identical and count>1, use '(X)n' form with sorted positions
        positions = [pos_labels[p] for p in pos_list]
        positions_sorted = sorted(positions, key=lambda x: int(x))
        if len(set(subs_list)) == 1 and len(subs_list) > 1:
            return f"{','.join(positions_sorted)}-({subs_list[0]}){len(subs_list)}-phen"
        else:
            # Mixed: e.g., '3,8-NO2,CHO-phen'
            return f"{','.join(positions_sorted)}-{'/'.join(subs_list)}-phen"

    # Modes enumeration
    mode_to_counts = {
        "mono": 1,
        "di": 2,
        "tri": 3,
        "tetra": 4,
    }

    for mode in modes:
        count = mode_to_counts[mode]
        # Position combinations
        pos_combos = list(itertools.combinations(sites, count))
        rng.shuffle(pos_combos)
        subs_list = list(sub_names)
        rng.shuffle(subs_list)
        for pos_idxs in pos_combos:
            # Substituent combinations
            if identical_only:
                subs_combos_iter = ([name] * count for name in subs_list)
            else:
                subs_product = itertools.product(subs_list, repeat=count)
                subs_combos_iter = subs_product
            for subs_combo in subs_combos_iter:
                # Build molecule
                mol = Chem.Mol(base)
                try:
                    for ring_idx, sub_name in zip(pos_idxs, subs_combo):
                        mol = attach_substituent(mol, ring_idx, sub_attach[sub_name])
                except Exception as e:
                    # Skip invalid combos (valence or sanitization issues)
                    continue
                smiles = Chem.MolToSmiles(mol, canonical=True)
                sys = system_label(list(pos_idxs), list(subs_combo))
                dv_c_sum = None
                # If all substituents map to DV_C (by name), sum them
                dv_values = []
                for s in subs_combo:
                    # Priority: refined mapping > e_int table > defaults
                    if s in dv_from_mapping:
                        dv_values.append(dv_from_mapping[s])
                    elif s in dv_from_eint:
                        dv_values.append(dv_from_eint[s])
                    elif s in DEFAULT_DV_C:
                        dv_values.append(DEFAULT_DV_C[s])
                    else:
                        dv_values.append(None)
                if all(v is not None for v in dv_values):
                    dv_c_sum = sum(v for v in dv_values if v is not None)
                rows.append({
                    "system": sys,
                    "smiles": smiles,
                    "sites": ",".join(pos_labels[i] for i in pos_idxs),
                    "substituents": ",".join(subs_combo),
                    "base_smiles": BASE_SMILES,
                    "n_subs": str(count),
                    "dv_c_sum": f"{dv_c_sum:.2f}" if dv_c_sum is not None else "",
                    "notes": "",
                })
                if max_rows is not None and len(rows) >= max_rows:
                    return rows
    return rows


def write_csv(path: str, rows: List[Dict[str, str]]):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fieldnames = ["system", "smiles", "sites", "substituents", "base_smiles", "n_subs", "dv_c_sum", "notes"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    args = parse_args()
    rows = generate_derivatives(
        args.subs,
        args.modes,
        args.eint,
        args.identical_only,
        args.max_rows,
        args.shuffle_seed,
    )
    write_csv(args.out, rows)
    print(f"Wrote {len(rows)} derivatives to {args.out}")


if __name__ == "__main__":
    main()
