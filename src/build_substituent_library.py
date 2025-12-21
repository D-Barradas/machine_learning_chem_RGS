import csv
import statistics
from collections import defaultdict, Counter
from pathlib import Path

DATA_PATHS = [
    Path("data/dv_values_by_cid.csv"),
    Path("data/dv_values.csv"),
    Path("data/substituent_selfies_propagated.csv"),
    Path("data/substituent_selfies_mapped.csv"),
    Path("data/substituent_selfies_mapped_pubchempy.csv"),
    Path("data/substituent_selfies_enriched.csv"),
    Path("data/substituent_selfies.csv"),
]
OUT_PATH = Path("data/substituent_dv_mappings.csv")


def normalize_name(name: str) -> str:
    if name is None:
        return ""
    # Trim, collapse whitespace, standardize simple variants
    n = " ".join(name.strip().split())
    # Remove trailing punctuation or commas often present in source
    n = n.strip(",.;")
    return n


def load_rows(limit: int | None = None):
    for path in DATA_PATHS:
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                if limit is not None and i >= limit:
                    break
                yield row


def aggregate():
    # Group per substituent
    dv_c_vals: dict[str, list[float]] = defaultdict(list)
    dv_n_vals: dict[str, list[float]] = defaultdict(list)
    smiles_counter: dict[str, Counter[str]] = defaultdict(Counter)
    count_counter: dict[str, int] = defaultdict(int)

    skipped = 0
    for row in load_rows():
        sub_raw = row.get("substituent", "")
        sub = normalize_name(sub_raw)
        if not sub:
            skipped += 1
            continue
        # Skip placeholder H (unsubstituted) entries
        if sub.upper() in {"H", "HYDROXIDE"}:
            continue

        dv_c_s = row.get("DV_C", "")
        dv_n_s = row.get("DV_n", "")
        smi = row.get("canonical_smiles", "").strip()

        try:
            dv_c = float(dv_c_s)
            dv_n = float(dv_n_s)
        except Exception:
            skipped += 1
            continue

        dv_c_vals[sub].append(dv_c)
        dv_n_vals[sub].append(dv_n)
        if smi:
            smiles_counter[sub][smi] += 1
        count_counter[sub] += 1

    # Build aggregated records
    records = []
    for sub, cvals in dv_c_vals.items():
        nvals = dv_n_vals.get(sub, [])
        if not cvals or not nvals:
            continue
        dv_c_med = statistics.median(cvals)
        dv_n_med = statistics.median(nvals)
        smi = ""
        if smiles_counter[sub]:
            # pick most common canonical_smiles
            smi = smiles_counter[sub].most_common(1)[0][0]
        records.append(
            {
                "substituent": sub,
                "canonical_smiles": smi,
                "DV_C": round(dv_c_med, 3),
                "DV_n": round(dv_n_med, 3),
                "examples": count_counter[sub],
            }
        )

    # Sort by examples desc, then by substituent name
    records.sort(key=lambda r: (-r["examples"], r["substituent"]))
    return records, skipped


def write_out(records):
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "substituent",
            "canonical_smiles",
            "DV_C",
            "DV_n",
            "examples",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def main():
    if not any(p.exists() for p in DATA_PATHS):
        raise FileNotFoundError(
            f"None of the sources found: {', '.join(str(p) for p in DATA_PATHS)}"
        )
    records, skipped = aggregate()
    write_out(records)
    print(
        f"Wrote {len(records)} mappings to {OUT_PATH} (skipped {skipped} rows)"
    )


if __name__ == "__main__":
    main()
