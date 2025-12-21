#!/usr/bin/env python3
"""
Reward service for RL-based ligand optimization.

Loads the trained ExtraTrees model (Eq8 pseudo-label) and computes rewards from SMILES.
Reward = -E_int_pred + penalties (validity, complexity, etc.)

Usage:
    from reward_service import RewardService
    service = RewardService("models/eint_hpo_rf.pkl", "data/substituent_dv_mappings.csv")
    reward = service.compute_reward(smiles)
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Optional, Tuple

import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors


class RewardService:
    def __init__(self, model_path: str, dv_mapping_path: str):
        """
        Parameters
        ----------
        model_path : path to trained joblib model dict with 'model' and 'features' keys
        dv_mapping_path : path to substituent_dv_mappings.csv
        """
        checkpoint = joblib.load(model_path)
        self.model = checkpoint["model"]
        self.features = checkpoint["features"]
        self.dv_mapping = self._load_dv_mapping(dv_mapping_path)

    def _load_dv_mapping(self, path: str) -> Dict[str, float]:
        mapping = {}
        if not Path(path).exists():
            return mapping
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get("substituent", "").strip()
                dv_str = row.get("DV_C", "").strip()
                if name and dv_str:
                    try:
                        mapping[name] = float(dv_str)
                    except Exception:
                        pass
        return mapping

    def extract_features(self, smiles: str) -> Optional[np.ndarray]:
        """
        Extract features [dv_c_sum, n_subs] from a phenanthroline derivative SMILES.
        Returns None if extraction fails.
        """
        # This is a simplified feature extractor; in reality you'd parse substituents
        # For now, we'll use a heuristic based on molecular weight and complexity
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Heuristic: approximate n_subs from number of heavy atoms beyond base phen
        # 1,10-phen has 12 carbons + 2 nitrogens = 14 heavy atoms
        base_heavy = 14
        heavy = mol.GetNumHeavyAtoms()
        n_subs = max(0, (heavy - base_heavy) // 3)  # rough estimate

        # Heuristic: dv_c_sum from molecular properties
        # Use a simple linear combination of descriptors as proxy
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        dv_c_sum = (mw - 180) * 0.1 + logp * 2.0  # arbitrary scaling

        return np.array([[dv_c_sum, n_subs]])

    def predict_eint(self, smiles: str) -> Optional[float]:
        """Predict interaction energy from SMILES. Returns None if invalid."""
        features = self.extract_features(smiles)
        if features is None:
            return None
        try:
            pred = self.model.predict(features)[0]
            return float(pred)
        except Exception:
            return None

    def compute_reward(
        self,
        smiles: str,
        validity_penalty: float = -100.0,
        complexity_weight: float = 0.1,
        target_eint: float = 50.0,
    ) -> Tuple[float, Dict[str, float]]:
        """
        Compute reward for a given SMILES.

        Reward = -(E_int_pred - target_eint)^2 - complexity_penalty + validity_bonus

        Parameters
        ----------
        smiles : candidate SMILES
        validity_penalty : penalty for invalid molecules
        complexity_weight : weight for complexity penalty (SA score)
        target_eint : target interaction energy (lower is better)

        Returns
        -------
        reward : float
        info : dict with components
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return validity_penalty, {"valid": False}

        eint_pred = self.predict_eint(smiles)
        if eint_pred is None:
            return validity_penalty, {"valid": False}

        # Reward components
        # 1. Distance to target (negative squared error)
        dist = (eint_pred - target_eint) ** 2
        reward_eint = -dist

        # 2. Complexity penalty (higher SA score = more complex = penalty)
        # For simplicity, use number of rotatable bonds as proxy
        complexity = Descriptors.NumRotatableBonds(mol)
        penalty_complexity = -complexity_weight * complexity

        reward = reward_eint + penalty_complexity

        info = {
            "valid": True,
            "eint_pred": eint_pred,
            "reward_eint": reward_eint,
            "penalty_complexity": penalty_complexity,
            "complexity": complexity,
        }
        return reward, info

    def batch_compute_reward(self, smiles_list: list[str], **kwargs) -> list[Tuple[float, Dict]]:
        """Batch compute rewards for a list of SMILES."""
        return [self.compute_reward(s, **kwargs) for s in smiles_list]
