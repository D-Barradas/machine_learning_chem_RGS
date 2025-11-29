#!/usr/bin/env python3
"""
PoC SELFIES tokenizer + tiny SELFIES-aware regressor (CLI)

Usage example:
  python src/selfies_tokenizer_poc.py \
    --train data/hf_dataset_small_train.csv \
    --val data/hf_dataset_small_val.csv \
    --test data/hf_dataset_small_test.csv \
    --target DV_C \
    --out models/poc_selfies/poc_selfies.pth

This script is a minimal, self-contained proof-of-concept. It does not aim to be production-ready.
"""

import argparse
import re
import os
import math
import json
from collections import Counter

import pandas as pd
import numpy as np

import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader

from sklearn.metrics import mean_squared_error

# Try to import selfies if available
try:
    import selfies as sf
    HAVE_SELFIES = True
except Exception:
    HAVE_SELFIES = False

BRACKET_TOKEN_RE = re.compile(r"\[[^\]]+\]")


def split_selfies(s):
    if not isinstance(s, str):
        return []
    # prefer library split if available
    if HAVE_SELFIES and hasattr(sf, "split_selfies"):
        try:
            return sf.split_selfies(s)
        except Exception:
            pass
    return BRACKET_TOKEN_RE.findall(s)


def build_vocab(series_of_selfies, min_freq=1, specials=None):
    if specials is None:
        specials = ["<PAD>", "<UNK>"]
    freq = Counter()
    for s in series_of_selfies.dropna():
        freq.update(split_selfies(s))
    tokens = [t for t, c in freq.items() if c >= min_freq]
    itos = specials + sorted(tokens)
    stoi = {t: i for i, t in enumerate(itos)}
    return stoi, itos


def encode_selfies(s, stoi, max_len=128):
    PAD_ID = stoi.get("<PAD>")
    UNK_ID = stoi.get("<UNK>")
    tokens = split_selfies(s)
    ids = [stoi.get(t, UNK_ID) for t in tokens]
    if len(ids) >= max_len:
        ids = ids[:max_len]
    pad_len = max_len - len(ids)
    return [PAD_ID] * pad_len + ids


class SelfiesDataset(Dataset):
    def __init__(self, df, stoi, max_len=128, target_col="DV_C"):
        self.texts = df["text"].fillna("").tolist()
        if target_col in df.columns:
            self.targets = df[target_col].astype(float).tolist()
        else:
            self.targets = [0.0] * len(self.texts)
        self.stoi = stoi
        self.max_len = max_len

    def __len__(self):
        return len(self.texts)

    def __getitem__(self, idx):
        ids = encode_selfies(self.texts[idx], self.stoi, max_len=self.max_len)
        ids = torch.tensor(ids, dtype=torch.long)
        target = torch.tensor(self.targets[idx], dtype=torch.float)
        return ids, target


class TinySelfiesRegressor(nn.Module):
    def __init__(self, vocab_size, embed_dim=64, pad_id=0):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, embed_dim, padding_idx=pad_id)
        self.linear = nn.Linear(embed_dim, 1)

    def forward(self, input_ids):
        embeds = self.embed(input_ids)  # (b, s, e)
        mask = (input_ids != self.embed.padding_idx).unsqueeze(-1).float()
        summed = (embeds * mask).sum(1)
        denom = mask.sum(1).clamp(min=1e-6)
        mean_pooled = summed / denom
        out = self.linear(mean_pooled).squeeze(-1)
        return out


def run_epoch(model, loader, loss_fn, optimizer=None, device="cpu"):
    training = optimizer is not None
    model.train() if training else model.eval()
    losses = []
    preds_all = []
    labels_all = []
    for xb, yb in loader:
        xb = xb.to(device)
        yb = yb.to(device)
        out = model(xb)
        loss = loss_fn(out, yb)
        if training:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        losses.append(loss.item())
        preds_all.extend(out.detach().cpu().numpy().tolist())
        labels_all.extend(yb.detach().cpu().numpy().tolist())
    mse = mean_squared_error(labels_all, preds_all) if len(preds_all) > 0 else float("nan")
    return np.mean(losses) if losses else float("nan"), math.sqrt(mse)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--train", required=True)
    p.add_argument("--val", required=True)
    p.add_argument("--test", required=True)
    p.add_argument("--target", default="DV_C")
    p.add_argument("--max_len", type=int, default=64)
    p.add_argument("--embed_dim", type=int, default=64)
    p.add_argument("--batch_size", type=int, default=8)
    p.add_argument("--epochs", type=int, default=5)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--out", default="models/poc_selfies/poc_selfies.pth")
    p.add_argument("--min_freq", type=int, default=1)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()
    torch.manual_seed(args.seed)

    train_df = pd.read_csv(args.train)
    val_df = pd.read_csv(args.val)
    test_df = pd.read_csv(args.test)

    # Build vocab from combined splits
    all_texts = pd.concat([train_df["text"], val_df["text"], test_df["text"]])
    stoi, itos = build_vocab(all_texts, min_freq=args.min_freq)
    PAD_ID = stoi.get("<PAD>")

    # Datasets
    train_ds = SelfiesDataset(train_df, stoi, max_len=args.max_len, target_col=args.target)
    val_ds = SelfiesDataset(val_df, stoi, max_len=args.max_len, target_col=args.target)
    test_ds = SelfiesDataset(test_df, stoi, max_len=args.max_len, target_col=args.target)

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size)
    test_loader = DataLoader(test_ds, batch_size=args.batch_size)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = TinySelfiesRegressor(vocab_size=len(itos), embed_dim=args.embed_dim, pad_id=PAD_ID)
    model.to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    loss_fn = nn.MSELoss()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    print(f"Vocab size: {len(itos)} | train/val/test: {len(train_ds)}/{len(val_ds)}/{len(test_ds)} | device: {device}")

    for epoch in range(args.epochs):
        train_loss, train_rmse = run_epoch(model, train_loader, loss_fn, optimizer=opt, device=device)
        val_loss, val_rmse = run_epoch(model, val_loader, loss_fn, optimizer=None, device=device)
        print(f"Epoch {epoch+1}/{args.epochs} - train_loss={train_loss:.4f} rmse={train_rmse:.4f} | val_rmse={val_rmse:.4f}")

    test_loss, test_rmse = run_epoch(model, test_loader, loss_fn, optimizer=None, device=device)
    print("Test RMSE:", test_rmse)

    # Save model and vocab
    out_dict = {
        "model_state_dict": model.state_dict(),
        "stoi": stoi,
        "itos": itos,
        "args": vars(args),
    }
    torch.save(out_dict, args.out)
    print("Saved model to", args.out)


if __name__ == "__main__":
    main()
