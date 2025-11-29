#!/usr/bin/env python3
"""Build a HuggingFace-style dataset CSV from SELFIES and DV labels.

Creates `text` (SELFIES), `DV_C`, `DV_n` columns and train/val/test splits.

Usage examples:
  # small mode from propagated small file
  python src/build_hf_dataset.py --input data/substituent_selfies_propagated.csv --output-prefix data/hf_dataset_small --mode small

  # large mode from dv_values_by_cid.csv (will convert SMILES -> SELFIES; may take long)
  python src/build_hf_dataset.py --input data/dv_values_by_cid.csv --output-prefix data/hf_dataset_large --mode large
"""
import argparse
from pathlib import Path
import sys
import pandas as pd


def to_selfies(smiles):
    try:
        import selfies
    except Exception:
        raise RuntimeError('selfies library not available; install with `pip install selfies`')
    try:
        return selfies.encoder(smiles) if isinstance(smiles, str) and smiles.strip() else None
    except Exception:
        return None


def build_from_small(inp_path: Path, out_prefix: Path, test_size=0.15, val_size=0.15, seed=42):
    df = pd.read_csv(inp_path)
    # ensure selfies column
    if 'selfies' not in df.columns or df['selfies'].isna().all():
        # try to compute from canonical_smiles
        if 'canonical_smiles' in df.columns:
            df['selfies'] = df['canonical_smiles'].apply(lambda s: to_selfies(s) if pd.notna(s) else None)
        else:
            raise RuntimeError('No `selfies` or `canonical_smiles` column present in small input')

    # keep rows where selfies exists and at least one DV label present
    df['DV_C'] = pd.to_numeric(df['DV_C'], errors='coerce') if 'DV_C' in df.columns else pd.NA
    df['DV_n'] = pd.to_numeric(df['DV_n'], errors='coerce') if 'DV_n' in df.columns else pd.NA
    mask = df['selfies'].notna() & (df['DV_C'].notna() | df['DV_n'].notna())
    df = df.loc[mask].copy()
    if df.empty:
        raise RuntimeError('No rows with SELFIES and DV labels found in input')

    out_all = out_prefix.with_suffix('.csv')
    # write full dataset
    hf = df[['selfies','DV_C','DV_n']].rename(columns={'selfies':'text'})
    hf.to_csv(out_all, index=False)

    # produce splits
    from sklearn.model_selection import train_test_split

    train_val, test = train_test_split(hf, test_size=test_size, random_state=seed)
    val_relative = val_size / (1.0 - test_size)
    train, val = train_test_split(train_val, test_size=val_relative, random_state=seed)

    train.to_csv(str(out_prefix) + '_train.csv', index=False)
    val.to_csv(str(out_prefix) + '_val.csv', index=False)
    test.to_csv(str(out_prefix) + '_test.csv', index=False)

    print('Wrote:', out_all)
    print('Splits:', str(out_prefix) + '_train.csv', str(out_prefix) + '_val.csv', str(out_prefix) + '_test.csv')
    print('Counts -> all:', len(hf), 'train:', len(train), 'val:', len(val), 'test:', len(test))


def build_from_large(inp_path: Path, out_prefix: Path, test_size=0.15, val_size=0.15, seed=42):
    df = pd.read_csv(inp_path)
    if 'canonical_smiles' not in df.columns:
        raise RuntimeError('Large mode requires `canonical_smiles` column in input')

    # convert SMILES -> SELFIES for unique SMILES to speed up
    unique_smiles = df['canonical_smiles'].dropna().astype(str).str.strip().unique()
    print('Unique SMILES to convert:', len(unique_smiles))
    s_map = {}
    for i, s in enumerate(unique_smiles, start=1):
        s_map[s] = to_selfies(s)
        if i % 1000 == 0:
            print('Converted', i, 'SMILES')

    df['selfies'] = df['canonical_smiles'].map(s_map)

    df['DV_C'] = pd.to_numeric(df['DV_C'], errors='coerce') if 'DV_C' in df.columns else pd.NA
    df['DV_n'] = pd.to_numeric(df['DV_n'], errors='coerce') if 'DV_n' in df.columns else pd.NA
    mask = df['selfies'].notna() & (df['DV_C'].notna() | df['DV_n'].notna())
    hf = df.loc[mask, ['selfies','DV_C','DV_n']].rename(columns={'selfies':'text'})
    out_all = out_prefix.with_suffix('.csv')
    hf.to_csv(out_all, index=False)

    # splits
    from sklearn.model_selection import train_test_split
    train_val, test = train_test_split(hf, test_size=test_size, random_state=seed)
    val_relative = val_size / (1.0 - test_size)
    train, val = train_test_split(train_val, test_size=val_relative, random_state=seed)

    train.to_csv(str(out_prefix) + '_train.csv', index=False)
    val.to_csv(str(out_prefix) + '_val.csv', index=False)
    test.to_csv(str(out_prefix) + '_test.csv', index=False)

    print('Wrote large dataset and splits. Counts -> all:', len(hf))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, help='Input CSV (small propagated or dv_values_by_cid.csv)')
    ap.add_argument('--output-prefix', required=True, help='Output prefix for files, e.g., data/hf_dataset_small')
    ap.add_argument('--mode', choices=['small','large'], default='small', help='small: use existing selfies; large: convert SMILES->SELFIES')
    ap.add_argument('--test-size', type=float, default=0.15)
    ap.add_argument('--val-size', type=float, default=0.15)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    inp = Path(args.input)
    outp = Path(args.output_prefix)
    if not inp.exists():
        print('Input not found:', inp)
        raise SystemExit(1)

    if args.mode == 'small':
        build_from_small(inp, outp, args.test_size, args.val_size, args.seed)
    else:
        build_from_large(inp, outp, args.test_size, args.val_size, args.seed)


if __name__ == '__main__':
    main()
