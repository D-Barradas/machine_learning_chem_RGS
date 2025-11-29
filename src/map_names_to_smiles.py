#!/usr/bin/env python3
"""Map substituent names to canonical SMILES using PubChem PUG-REST.

Usage:
  python src/map_names_to_smiles.py --input data/substituent_selfies_enriched.csv --output data/substituent_selfies_mapped.csv

For rows missing `canonical_smiles`, the script will query PubChem by `substituent` (or `Name`) and fill `canonical_smiles` when found.
If the `selfies` library is available and `--write-selfies` is passed, the script will also encode SELFIES for newly found SMILES.
"""
import argparse
from pathlib import Path
import time
import sys
import pandas as pd

try:
    from urllib.parse import quote_plus
    import requests
except Exception as e:
    print('Missing standard library modules?', e, file=sys.stderr)
    raise

# Prefer pubchempy when available for robust name/inchikey/synonym lookups
try:
    import pubchempy as pcp
except Exception:
    pcp = None


def get_smiles_from_pubchem(name, timeout=10):
    """Resolve a name (or InChIKey) to canonical SMILES using pubchempy when available,
    otherwise fall back to PUG-REST. Returns None if not found.
    """
    if not isinstance(name, str) or not name.strip():
        return None
    # detect possible InChIKey (27 chars with hyphen parts)
    s = name.strip()
    if pcp is not None:
        try:
            # try by InChIKey first
            if len(s) >= 27 and '-' in s:
                comps = pcp.get_compounds(s, 'inchikey')
                if comps:
                    return comps[0].canonical_smiles
            # try by CID if numeric
            if s.isdigit():
                comps = pcp.get_compounds(int(s), 'cid')
                if comps:
                    return comps[0].canonical_smiles
            # try by name
            comps = pcp.get_compounds(s, 'name')
            if comps:
                return comps[0].canonical_smiles
            # try synonyms search (returns list)
            # pubchempy doesn't have a direct 'synonym' query return type for get_compounds,
            # but searching by name above often handles common synonyms.
        except Exception:
            pass

    # fallback to PUG-REST
    try:
        q = quote_plus(s)
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/CanonicalSMILES/JSON'
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200:
            j = r.json()
            props = j.get('PropertyTable', {}).get('Properties', [])
            if props:
                smi = props[0].get('CanonicalSMILES')
                return smi
    except Exception:
        return None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, help='Input CSV (enriched)')
    ap.add_argument('--output', required=True, help='Output CSV with mapped SMILES')
    ap.add_argument('--write-selfies', action='store_true', help='Encode SELFIES when possible')
    ap.add_argument('--pause', type=float, default=0.2, help='Pause between PubChem requests (s)')
    args = ap.parse_args()

    inp = Path(args.input)
    outp = Path(args.output)
    if not inp.exists():
        print(f'Input {inp} not found', file=sys.stderr)
        raise SystemExit(1)

    df = pd.read_csv(inp)
    df_columns = list(df.columns)
    # define candidate name columns
    name_cols = [c for c in ['substituent', 'Name', 'name'] if c in df.columns]

    missing_idx = df[df['canonical_smiles'].isna() | (df['canonical_smiles'].astype(str).str.strip()=='')].index
    print(f'Rows missing canonical_smiles: {len(missing_idx)}')

    use_selfies = False
    if args.write_selfies:
        try:
            import selfies
            use_selfies = True
        except Exception:
            print('`selfies` not available; continuing without encoding SELFIES', file=sys.stderr)
            use_selfies = False

    updated = 0
    for i in missing_idx:
        # pick a name to query
        name = None
        for nc in name_cols:
            val = df.at[i, nc]
            if isinstance(val, str) and val.strip():
                name = val.strip()
                break
        if not name:
            continue
        smi = get_smiles_from_pubchem(name)
        time.sleep(args.pause)
        if smi:
            df.at[i, 'canonical_smiles'] = smi
            if use_selfies:
                try:
                    df.at[i, 'selfies'] = selfies.encoder(smi)
                except Exception:
                    df.at[i, 'selfies'] = df.at[i, 'selfies']
            updated += 1
            if updated % 10 == 0:
                print(f'Updated {updated} rows so far...')

    df.to_csv(outp, index=False)
    print(f'Wrote {outp} with {len(df)} rows; updated {updated} SMILES')


if __name__ == '__main__':
    main()
