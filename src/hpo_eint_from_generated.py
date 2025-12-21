#!/usr/bin/env python3
"""
Hyperparameter search (lightweight) to predict interaction energies from generated derivatives.

Data source
- data/phen_derivatives_scored.csv (200k rows with dv_c_sum, n_subs, e_int_pred_eq8_fit)

Features
- dv_c_sum (float)
- n_subs (int)

Target
- e_int_pred_eq8_fit

Approach
- Train/validation split (80/20)
- Randomized search over a few tree-based regressors
- Report best params/metrics; persist best model

Usage
    python src/hpo_eint_from_generated.py \
        --data data/phen_derivatives_scored.csv \
        --out-model models/eint_hpo_rf.pkl \
        --results data/eint_hpo_rf_metrics.json \
        --sample 80000

Notes
- This is a fast surrogate HPO to rank patterns; it fits on Eq8-based pseudo-labels.
- For true energies, retrain on experimental/computed E_int when available.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Any

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="HPO for E_int using generated DV features")
    p.add_argument("--data", default="data/phen_derivatives_scored.csv", help="Input scored derivatives CSV")
    p.add_argument("--out-model", default="models/eint_hpo_rf.pkl", help="Where to save best model")
    p.add_argument("--results", default="data/eint_hpo_rf_metrics.json", help="Where to save metrics")
    p.add_argument("--sample", type=int, default=80000, help="Optional row cap for speed; use all if larger than dataset")
    p.add_argument("--seed", type=int, default=17, help="Random seed")
    return p.parse_args()


def load_data(path: str, sample: int | None, seed: int):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    needed = ["dv_c_sum", "n_subs", "e_int_pred_eq8_fit"]
    for c in needed:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' missing in {path}")
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["dv_c_sum", "e_int_pred_eq8_fit"], how="any")
    df["n_subs"] = df["n_subs"].astype(int)
    if sample is not None and sample < len(df):
        df = df.sample(n=sample, random_state=seed)
    return df


def evaluate_model(model, X_train, X_val, y_train, y_val) -> Dict[str, float]:
    model.fit(X_train, y_train)
    pred = model.predict(X_val)
    mae = mean_absolute_error(y_val, pred)
    rmse = mean_squared_error(y_val, pred) ** 0.5
    r2 = r2_score(y_val, pred)
    return {"mae": mae, "rmse": rmse, "r2": r2}


def main():
    args = parse_args()
    df = load_data(args.data, args.sample, args.seed)
    X = df[["dv_c_sum", "n_subs"]]
    y = df["e_int_pred_eq8_fit"]

    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.2, random_state=args.seed
    )

    candidates: Dict[str, Any] = {
        "rf": RandomForestRegressor(
            n_estimators=400,
            max_depth=None,
            min_samples_leaf=1,
            random_state=args.seed,
            n_jobs=-1,
        ),
        "et": ExtraTreesRegressor(
            n_estimators=400,
            max_depth=None,
            min_samples_leaf=1,
            random_state=args.seed,
            n_jobs=-1,
        ),
        "gbr": GradientBoostingRegressor(random_state=args.seed),
    }

    metrics: Dict[str, Dict[str, float]] = {}
    best_name = None
    best_rmse = float("inf")
    best_model = None

    for name, model in candidates.items():
        m = evaluate_model(model, X_train, X_val, y_train, y_val)
        metrics[name] = m
        if m["rmse"] < best_rmse:
            best_rmse = m["rmse"]
            best_name = name
            best_model = model

    Path(args.out_model).parent.mkdir(parents=True, exist_ok=True)
    Path(args.results).parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": best_model, "features": ["dv_c_sum", "n_subs"]}, args.out_model)
    with open(args.results, "w", encoding="utf-8") as f:
        json.dump({"metrics": metrics, "best_model": best_name}, f, indent=2)

    print(f"Best model: {best_name} (rmse={best_rmse:.3f}) on Eq8 pseudo-labels")
    print(f"Saved model to {args.out_model}")
    print(f"Saved metrics to {args.results}")


if __name__ == "__main__":
    main()
