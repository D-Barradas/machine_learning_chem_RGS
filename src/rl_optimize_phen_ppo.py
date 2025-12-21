#!/usr/bin/env python3
"""
PPO-based RL optimization for phenanthroline derivatives using SELFIES.

Trains a policy network to generate SELFIES sequences that maximize predicted binding affinity.

Architecture:
- Policy: LSTM over SELFIES tokens
- Reward: ExtraTrees predictor (Eq8 pseudo-label)
- Algorithm: PPO with KL regularization

Usage:
    python src/rl_optimize_phen_ppo.py \
        --model models/eint_hpo_rf.pkl \
        --dv-mapping data/substituent_dv_mappings.csv \
        --smiles-data data/phen_derivatives_scored.csv \
        --out-dir results/rl_ppo \
        --episodes 200 \
        --batch-size 32

Outputs:
- Trained policy checkpoint
- CSV of top-k candidates with SMILES/SELFIES/reward
- Training logs (rewards per episode)
"""
from __future__ import annotations

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Categorical

from reward_service import RewardService
from selfies_tokenizer import SELFIESTokenizer


class PolicyNetwork(nn.Module):
    """LSTM-based policy for SELFIES token generation."""

    def __init__(self, vocab_size: int, embedding_dim: int = 128, hidden_dim: int = 256):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.fc = nn.Linear(hidden_dim, vocab_size)

    def forward(self, x: torch.Tensor, hidden=None) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        x: (batch, seq_len)
        Returns: logits (batch, seq_len, vocab_size), hidden
        """
        emb = self.embedding(x)
        out, hidden = self.lstm(emb, hidden)
        logits = self.fc(out)
        return logits, hidden

    def sample(self, tokenizer: SELFIESTokenizer, max_len: int = 50, temperature: float = 1.0) -> Tuple[List[int], float]:
        """Sample a sequence from the policy."""
        device = next(self.parameters()).device
        tokens = [tokenizer.sos_id]
        log_probs = []
        hidden = None

        for _ in range(max_len):
            x = torch.tensor([[tokens[-1]]], device=device)
            logits, hidden = self.forward(x, hidden)
            logits = logits[:, -1, :] / temperature
            probs = F.softmax(logits, dim=-1)
            dist = Categorical(probs)
            action = dist.sample()
            log_prob = dist.log_prob(action)

            tokens.append(action.item())
            log_probs.append(log_prob.item())

            if tokens[-1] == tokenizer.eos_id:
                break

        return tokens, sum(log_probs)


class PPOTrainer:
    """PPO trainer for molecule generation."""

    def __init__(
        self,
        policy: PolicyNetwork,
        tokenizer: SELFIESTokenizer,
        reward_service: RewardService,
        lr: float = 1e-4,
        gamma: float = 0.99,
        epsilon: float = 0.2,
        entropy_coef: float = 0.01,
    ):
        self.policy = policy
        self.tokenizer = tokenizer
        self.reward_service = reward_service
        self.optimizer = torch.optim.Adam(policy.parameters(), lr=lr)
        self.gamma = gamma
        self.epsilon = epsilon
        self.entropy_coef = entropy_coef

    def collect_trajectories(self, batch_size: int, max_len: int = 50) -> List[Dict]:
        """Sample trajectories from current policy."""
        trajectories = []
        for _ in range(batch_size):
            tokens, log_prob_sum = self.policy.sample(self.tokenizer, max_len=max_len)
            selfies_str = self.tokenizer.decode_tokens(tokens)
            smiles = self.tokenizer.selfies_to_smiles(selfies_str)

            if smiles is None:
                reward = -100.0
                info = {"valid": False}
            else:
                reward, info = self.reward_service.compute_reward(smiles)

            trajectories.append({
                "tokens": tokens,
                "log_prob": log_prob_sum,
                "reward": reward,
                "smiles": smiles,
                "selfies": selfies_str,
                "info": info,
            })
        return trajectories

    def update_policy(self, trajectories: List[Dict], epochs: int = 4) -> Dict[str, float]:
        """PPO update step."""
        # Prepare data
        all_tokens = []
        old_log_probs = []
        rewards = []
        for traj in trajectories:
            all_tokens.append(traj["tokens"])
            old_log_probs.append(traj["log_prob"])
            rewards.append(traj["reward"])

        # Normalize rewards
        rewards = np.array(rewards)
        if rewards.std() > 1e-6:
            rewards = (rewards - rewards.mean()) / (rewards.std() + 1e-8)

        device = next(self.policy.parameters()).device
        losses = []

        for _ in range(epochs):
            for i, tokens in enumerate(all_tokens):
                if len(tokens) < 2:
                    continue

                # Compute new log probs
                x = torch.tensor([tokens[:-1]], device=device, dtype=torch.long)
                logits, _ = self.policy.forward(x)
                log_probs = F.log_softmax(logits, dim=-1)

                # Gather log probs for taken actions
                actions = torch.tensor([tokens[1:]], device=device, dtype=torch.long)
                selected_log_probs = log_probs.gather(2, actions.unsqueeze(-1)).squeeze(-1)
                new_log_prob = selected_log_probs.sum()

                # PPO ratio
                old_lp = old_log_probs[i]
                ratio = torch.exp(new_log_prob - old_lp)
                advantage = torch.tensor(rewards[i], device=device)

                # PPO clipped objective
                surr1 = ratio * advantage
                surr2 = torch.clamp(ratio, 1.0 - self.epsilon, 1.0 + self.epsilon) * advantage
                policy_loss = -torch.min(surr1, surr2)

                # Entropy bonus
                probs = F.softmax(logits, dim=-1)
                entropy = -(probs * log_probs).sum(dim=-1).mean()
                loss = policy_loss - self.entropy_coef * entropy

                self.optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.policy.parameters(), 1.0)
                self.optimizer.step()

                losses.append(loss.item())

        return {"loss": np.mean(losses) if losses else 0.0}

    def train(self, episodes: int, batch_size: int, save_dir: Path):
        """Main training loop."""
        save_dir.mkdir(parents=True, exist_ok=True)
        log_path = save_dir / "training_log.json"
        best_path = save_dir / "best_policy.pt"
        top_mols_path = save_dir / "top_molecules.csv"

        logs = []
        best_reward = -float("inf")
        top_molecules = []

        for episode in range(episodes):
            trajectories = self.collect_trajectories(batch_size)
            update_info = self.update_policy(trajectories)

            # Track metrics
            rewards = [t["reward"] for t in trajectories]
            valid_count = sum(1 for t in trajectories if t["info"].get("valid", False))
            mean_reward = np.mean(rewards)
            max_reward = np.max(rewards)

            logs.append({
                "episode": episode,
                "mean_reward": mean_reward,
                "max_reward": max_reward,
                "valid_ratio": valid_count / batch_size,
                "loss": update_info["loss"],
            })

            # Update top molecules
            for traj in trajectories:
                if traj["smiles"] is not None:
                    top_molecules.append({
                        "smiles": traj["smiles"],
                        "selfies": traj["selfies"],
                        "reward": traj["reward"],
                        "eint_pred": traj["info"].get("eint_pred", None),
                    })

            # Keep top 1000
            top_molecules = sorted(top_molecules, key=lambda x: x["reward"], reverse=True)[:1000]

            # Save best policy
            if max_reward > best_reward:
                best_reward = max_reward
                torch.save(self.policy.state_dict(), best_path)

            # Print progress
            if (episode + 1) % 10 == 0:
                print(f"Episode {episode+1}/{episodes} | Mean Reward: {mean_reward:.3f} | "
                      f"Max Reward: {max_reward:.3f} | Valid: {valid_count}/{batch_size}")

        # Save logs and top molecules
        with open(log_path, "w") as f:
            json.dump(logs, f, indent=2)

        with open(top_mols_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["smiles", "selfies", "reward", "eint_pred"])
            writer.writeheader()
            writer.writerows(top_molecules)

        print(f"Training complete. Best reward: {best_reward:.3f}")
        print(f"Logs saved to {log_path}")
        print(f"Top molecules saved to {top_mols_path}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="PPO-based RL optimization for phen derivatives")
    p.add_argument("--model", default="models/eint_hpo_rf.pkl", help="Trained reward model path")
    p.add_argument("--dv-mapping", default="data/substituent_dv_mappings.csv", help="DV_C mapping CSV")
    p.add_argument("--smiles-data", default="data/phen_derivatives_scored.csv", help="SMILES data for vocab")
    p.add_argument("--out-dir", default="results/rl_ppo", help="Output directory")
    p.add_argument("--episodes", type=int, default=200, help="Number of training episodes")
    p.add_argument("--batch-size", type=int, default=32, help="Trajectories per episode")
    p.add_argument("--max-len", type=int, default=50, help="Max SELFIES sequence length")
    p.add_argument("--lr", type=float, default=1e-4, help="Learning rate")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    return p.parse_args()


def main():
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)

    # Force CPU for stability (cuDNN issues on some systems)
    device = torch.device("cpu")
    print(f"Using device: {device}")

    # Initialize components
    print("Building tokenizer...")
    tokenizer = SELFIESTokenizer.from_smiles_file(args.smiles_data)
    print(f"Vocab size: {tokenizer.vocab_size()}")

    print("Loading reward service...")
    reward_service = RewardService(args.model, args.dv_mapping)

    print("Initializing policy network...")
    policy = PolicyNetwork(tokenizer.vocab_size()).to(device)

    print("Creating PPO trainer...")
    trainer = PPOTrainer(policy, tokenizer, reward_service, lr=args.lr)

    print(f"Starting training for {args.episodes} episodes...")
    trainer.train(args.episodes, args.batch_size, Path(args.out_dir))


if __name__ == "__main__":
    main()
