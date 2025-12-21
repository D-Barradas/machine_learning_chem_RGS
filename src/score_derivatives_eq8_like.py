#!/usr/bin/env python3
"""
Fit simple linear models E_int ≈ a + b * DV_C using Table3_eq8 rows,
then score generated derivatives using summed DV_C of substituents.

Inputs:
- data/e_int_consolidated.csv (provided)
- data/phen_derivatives_generated_fix.csv (from generator)

Outputs:
- data/phen_derivatives_scored.csv with added column e_int_pred_eq8_fit

Note: This approximates Eq.7/8 behavior for mixed substituents by using DV_C sum.
"""
import csv
import os
from typing import Dict, List, Tuple

import numpy as np
from sklearn.linear_model import LinearRegression

EINT_PATH = "data/e_int_consolidated.csv"
DERIV_PATH = "data/phen_derivatives_generated_fix.csv"
OUT_PATH = "data/phen_derivatives_scored.csv"


def read_table3(path: str) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
    xs: List[float] = []
    ys_mono: List[float] = []
    ys_di: List[float] = []
    ys_tetra: List[float] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("Source", "") != "Table3_eq8":
                continue
            try:
                dv = float(row.get("DV_C (kcal/mol)", "").strip())
            except Exception:
                continue
            xs.append(dv)
            # Collect y values if present
            def parse_float(key: str):
                v = row.get(key, "").strip()
                try:
                    return float(v)
                except Exception:
                    return None
            mono = parse_float("E_int_eq8_mono (kcal/mol)")
            di = parse_float("E_int_eq8_di (kcal/mol)")
            tetra = parse_float("E_int_eq8_tetra (kcal/mol)")
            ys_mono.append(mono if mono is not None else np.nan)
            ys_di.append(di if di is not None else np.nan)
            ys_tetra.append(tetra if tetra is not None else np.nan)
    x = np.array(xs).reshape(-1, 1)
    return x, {1: np.array(ys_mono), 2: np.array(ys_di), 4: np.array(ys_tetra)}


def fit_models(x: np.ndarray, y_map: Dict[int, np.ndarray]) -> Dict[int, LinearRegression]:
    models: Dict[int, LinearRegression] = {}
    for k, y in y_map.items():
        mask = ~np.isnan(y)
        xk = x[mask]
        yk = y[mask]
        if len(xk) < 3:
            continue
        m = LinearRegression().fit(xk, yk)
        models[k] = m
    return models


def read_dv_c_sum_from_derivatives(path: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(row)
    return rows


def write_scored(rows: List[Dict[str, str]], preds: List[float], out_path: str):
    fieldnames = list(rows[0].keys()) + ["e_int_pred_eq8_fit"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row, p in zip(rows, preds):
            row_out = dict(row)
            row_out["e_int_pred_eq8_fit"] = f"{p:.2f}" if p is not None else ""
            w.writerow(row_out)


def main():
    x, y_map = read_table3(EINT_PATH)
    models = fit_models(x, y_map)
    deriv_rows = read_dv_c_sum_from_derivatives(DERIV_PATH)
    preds: List[float] = []
    for row in deriv_rows:
        try:
            dv_sum = float((row.get("dv_c_sum") or "").strip())
        except Exception:
            dv_sum = None
        n_subs = int(row.get("n_subs", "1"))
        if dv_sum is None:
            preds.append(None)
            continue
        if n_subs in models:
            p = models[n_subs].predict(np.array([[dv_sum]]))[0]
            preds.append(float(p))
            continue
        # Fallback for tri-substitution (3): interpolate between di (2) and tetra (4)
        if n_subs == 3:
            p2 = models.get(2)
            p4 = models.get(4)
            if p2 and p4:
                p_di = float(p2.predict(np.array([[dv_sum]]))[0])
                p_tetra = float(p4.predict(np.array([[dv_sum]]))[0])
                preds.append((p_di + p_tetra) / 2.0)
            elif p2:
                preds.append(float(p2.predict(np.array([[dv_sum]]))[0]))
            elif p4:
                preds.append(float(p4.predict(np.array([[dv_sum]]))[0]))
            else:
                preds.append(None)
        else:
            preds.append(None)
    write_scored(deriv_rows, preds, OUT_PATH)
    print(f"Wrote scored derivatives to {OUT_PATH}")


if __name__ == "__main__":
    main()
