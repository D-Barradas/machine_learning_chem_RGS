#!/usr/bin/env python3
"""Prepare dataset by converting SMILES to SELFIES and collecting DV_C / DV_n labels.

Usage:
  python src/prepare_dataset_selfies.py --inputs data/preliminar_data_961.csv data/preliminar_data_0961.csv --output data/substituent_selfies.csv

This script tries to be robust to slight variations in column names (e.g., "ΔVC-m", "ΔVC -p").
"""
import argparse
from pathlib import Path
import sys
import pandas as pd


def find_dv_columns(columns):
    # Return candidates for DV_C (m) and DV_n (p) columns
    dv_m = None
    dv_p = None
    for c in columns:
        low = c.lower()
        # look for delta-V columns
        if 'Δvc' in c or 'Δvc' in low or 'dvc' in low or 'delta' in low:
            # prefer columns that mention m or p explicitly
            if 'm' in low and dv_m is None:
                dv_m = c
            elif 'p' in low and dv_p is None:
                dv_p = c
            else:
                if dv_m is None:
                    dv_m = c
    return dv_m, dv_p


def convert_files(input_files, output_path):
    try:
        import selfies
    except Exception as e:
        print('ERROR: the `selfies` library is required but not installed.', file=sys.stderr)
        print('Install with: pip install selfies', file=sys.stderr)
        raise

    rows = []
    for f in input_files:
        p = Path(f)
        if not p.exists():
            print(f'Warning: {p} does not exist — skipping', file=sys.stderr)
            continue
        print(f'Reading {p}...')
        df = pd.read_csv(p)
        dv_m_col, dv_p_col = find_dv_columns(df.columns)
        print('Detected DV columns:', dv_m_col, dv_p_col)

        for _, r in df.iterrows():
            substituent = r.get('substituent') or r.get('Substituent') or r.get('Name') or ''
            cid = r.get('cid') if 'cid' in df.columns else ''
            smiles = r.get('canonical_smiles') or r.get('smiles') or ''
            dv_m = r.get(dv_m_col) if dv_m_col in df.columns else None
            dv_p = r.get(dv_p_col) if (dv_p_col and dv_p_col in df.columns) else None

            selfies_str = None
            if isinstance(smiles, str) and smiles.strip():
                try:
                    selfies_str = selfies.encoder(smiles)
                except Exception:
                    selfies_str = None

            rows.append({
                'substituent': substituent,
                'cid': cid,
                'canonical_smiles': smiles,
                'selfies': selfies_str,
                'DV_C': dv_m,
                'DV_n': dv_p,
                'source_file': str(p)
            })

    out_df = pd.DataFrame(rows)
    out_path = Path(output_path)
    out_df.to_csv(out_path, index=False)
    print(f'Wrote {len(out_df)} rows to {out_path}')


def main():
    ap = argparse.ArgumentParser(description='Build substituent SELFIES dataset')
    ap.add_argument('--inputs', nargs='+', required=True, help='Input CSV files (preliminar_data... )')
    ap.add_argument('--output', required=True, help='Output CSV file, e.g., data/substituent_selfies.csv')
    args = ap.parse_args()
    convert_files(args.inputs, args.output)


if __name__ == '__main__':
    main()
