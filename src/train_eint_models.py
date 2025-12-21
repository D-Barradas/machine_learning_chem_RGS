#!/usr/bin/env python3
"""
Train and select a regression model to predict interaction energy (E_int) from DV descriptors.

Inputs
- data/e_int_consolidated.csv (default) with columns like:
  - DV_N (kcal/mol)
  - DV_C (kcal/mol)
  - E_int_calc (kcal/mol)  <-- target

Outputs
- JSON with CV metrics per model
- Serialized best model (joblib)

Usage
    python src/train_eint_models.py \
        --data data/e_int_consolidated.csv \
        --target "E_int_calc (kcal/mol)" \
        --results data/eint_model_scores.json \
        --out-model models/eint_best_model.pkl

Notes
- Drops rows missing the target and all features; imputes remaining missing feature values with median.
- Candidate models: linear, ridge, random forest, gradient boosting. Best picked by lowest mean RMSE.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Dict, List

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Train/selection loop for E_int prediction")
    p.add_argument("--data", default="data/e_int_consolidated.csv", help="Input CSV path")
    p.add_argument("--target", default="E_int_calc (kcal/mol)", help="Primary target column name")
    p.add_argument("--results", default="data/eint_model_scores.json", help="Where to write CV metrics")
    p.add_argument("--out-model", default="models/eint_best_model.pkl", help="Where to store best model")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--folds", type=int, default=5, help="KFold splits (reduces if dataset is small)")
    p.add_argument(
        "--models",
        nargs="+",
        default=["linear", "ridge", "rf", "gbr"],
        choices=["linear", "ridge", "rf", "gbr"],
        help="Which models to evaluate",
    )
    return p.parse_args()


def load_data(path: str, target_col: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Normalize column names for easier access
    df.columns = [c.strip() for c in df.columns]
    target_col = target_col.strip()
    if target_col not in df.columns:
        raise ValueError(f"Target column '{target_col}' not found in {path}")
    return df


def prepare_features(df: pd.DataFrame, target_col: str):
    feature_candidates = [
        "DV_N (kcal/mol)",
        "DV_C (kcal/mol)",
    ]
    features = [c for c in feature_candidates if c in df.columns]
    if not features:
        raise ValueError("No DV feature columns found in dataset")

    data = df[[c for c in features + [target_col]]].copy()
    for c in features + [target_col]:
        data[c] = pd.to_numeric(data[c], errors="coerce")

    # Drop rows without target and without all features
    data = data.dropna(subset=[target_col], how="any")
    # Drop feature columns that are entirely missing for this subset
    features_with_data = [c for c in features if data[c].notna().any()]
    if not features_with_data:
        raise ValueError(f"No usable feature columns for target '{target_col}' (all NaN)")
    data = data.dropna(subset=features_with_data, how="all")

    # If no rows remain, signal failure so caller can try alternative targets
    if len(data) == 0:
        raise ValueError(f"No rows available after filtering for target '{target_col}'")

    X = data[features_with_data]
    y = data[target_col]
    return X, y, features_with_data


def build_models(seed: int) -> Dict[str, Pipeline]:
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()

    models: Dict[str, Pipeline] = {
        "linear": Pipeline([("imputer", imputer), ("scaler", scaler), ("model", LinearRegression())]),
        "ridge": Pipeline([
            ("imputer", imputer),
            ("scaler", scaler),
            ("model", Ridge(alpha=1.0, random_state=seed)),
        ]),
        "rf": Pipeline([
            ("imputer", imputer),
            ("model", RandomForestRegressor(
                n_estimators=400, max_depth=None, random_state=seed, n_jobs=-1
            )),
        ]),
        "gbr": Pipeline([
            ("imputer", imputer),
            ("scaler", scaler),
            ("model", GradientBoostingRegressor(random_state=seed)),
        ]),
    }
    return models


def evaluate_model(name: str, model: Pipeline, X: pd.DataFrame, y: pd.Series, folds: int, seed: int):
    n_samples = len(X)
    k = min(folds, n_samples)
    if k < 2:
        # Fit and evaluate on training data (best-effort when data is extremely small)
        model.fit(X, y)
        y_pred = model.predict(X)
        mae = mean_absolute_error(y, y_pred)
        rmse = mean_squared_error(y, y_pred, squared=False)
        r2 = r2_score(y, y_pred)
        return {"mae": mae, "rmse": rmse, "r2": r2, "cv": False, "folds": k}
    cv = KFold(n_splits=k, shuffle=True, random_state=seed)
    y_pred = cross_val_predict(model, X, y, cv=cv, n_jobs=-1)
    mae = mean_absolute_error(y, y_pred)
    rmse = mean_squared_error(y, y_pred) ** 0.5
    r2 = r2_score(y, y_pred)
    return {"mae": mae, "rmse": rmse, "r2": r2, "cv": True, "folds": k}


def select_and_train(models: Dict[str, Pipeline], X: pd.DataFrame, y: pd.Series, folds: int, seed: int):
    metrics: Dict[str, Dict[str, float]] = {}
    for name, model in models.items():
        metrics[name] = evaluate_model(name, model, X, y, folds, seed)

    best = min(metrics.items(), key=lambda kv: kv[1]["rmse"])
    best_name = best[0]
    best_model = models[best_name]
    best_model.fit(X, y)
    return best_name, best_model, metrics


def main():
    args = parse_args()
    df = load_data(args.data, args.target)

    target_candidates: List[str] = [
        args.target,
        "E_int_pred_eq5 (kcal/mol)",
        "E_int_eq8_mono (kcal/mol)",
    ]

    X = y = features = None
    chosen_target = None
    last_error: Exception | None = None
    for tgt in target_candidates:
        try:
            X, y, features = prepare_features(df, tgt)
            chosen_target = tgt
            break
        except Exception as e:
            last_error = e
            continue

    if chosen_target is None or X is None or y is None:
        raise RuntimeError(f"Failed to prepare data for any target ({target_candidates}): {last_error}")
    models = build_models(args.seed)
    selected_models = {k: v for k, v in models.items() if k in args.models}

    best_name, best_model, metrics = select_and_train(selected_models, X, y, args.folds, args.seed)

    # Ensure output directories exist
    Path(args.out_model).parent.mkdir(parents=True, exist_ok=True)
    Path(args.results).parent.mkdir(parents=True, exist_ok=True)

    joblib.dump({"model": best_model, "features": features, "target": chosen_target}, args.out_model)
    with open(args.results, "w", encoding="utf-8") as f:
        json.dump({"metrics": metrics, "best_model": best_name, "features": features, "target": chosen_target}, f, indent=2)

    print(f"Target: {chosen_target}")
    print(f"Best model: {best_name} (rmse={metrics[best_name]['rmse']:.3f}, r2={metrics[best_name]['r2']:.3f})")
    print(f"Saved model to {args.out_model}")
    print(f"Saved metrics to {args.results}")


if __name__ == "__main__":
    main()
