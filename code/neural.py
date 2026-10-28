"""TARNet (Shalit et al. 2017) and DragonNet (Shi et al. 2019) in PyTorch, with
the frozen-plan (§4, E4/E5) specification:

* representation MLP 3 x 200 (ELU); two outcome heads 2 x 100 (ELU) -> linear;
  DragonNet adds a linear+sigmoid propensity head on the representation.
* inputs and outcome z-scored per cohort; predictions un-standardised.
* loss = factual MSE (+ 1.0 * BCE(propensity) for DragonNet); L2 1e-4 (weight decay).
* Adam lr 1e-3, batch 100, <= 500 epochs, early stopping on a seeded 30 %
  validation split (stratified on T): patience 30 epochs, min 50 epochs, best-val
  weights restored.
* K seeds; in-sample mu0(x), mu1(x) [, e(x)] for every unit.

catenets is not installed (JAX dependency -> install requires user approval), so
this is a faithful re-implementation with catenets SNet1/SNet2 layer sizes.
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

torch.set_num_threads(1)


def _mlp(n_in, n_hidden, n_layers, n_out=None):
    layers = []
    d = n_in
    for _ in range(n_layers):
        layers += [nn.Linear(d, n_hidden), nn.ELU()]
        d = n_hidden
    if n_out is not None:
        layers += [nn.Linear(d, n_out)]
    return nn.Sequential(*layers)


class TarNet(nn.Module):
    def __init__(self, p, n_units_r=200, n_layers_r=3, n_units_out=100, n_layers_out=2,
                 propensity_head=False):
        super().__init__()
        self.rep = _mlp(p, n_units_r, n_layers_r)
        self.head0 = _mlp(n_units_r, n_units_out, n_layers_out, 1)
        self.head1 = _mlp(n_units_r, n_units_out, n_layers_out, 1)
        self.prop = nn.Linear(n_units_r, 1) if propensity_head else None

    def forward(self, x):
        h = self.rep(x)
        y0 = self.head0(h).squeeze(-1)
        y1 = self.head1(h).squeeze(-1)
        logit = self.prop(h).squeeze(-1) if self.prop is not None else None
        return y0, y1, logit


def fit_predict(X, T, Y, seed: int, dragon: bool = False, lr: float = 1e-3, batch: int = 100,
                max_epochs: int = 500, patience: int = 30, min_epochs: int = 50,
                val_frac: float = 0.3, weight_decay: float = 1e-4, alpha_prop: float = 1.0):
    """Train one network; return (mu0, mu1, e_hat_or_None, info) for ALL units (in-sample)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    n, p = X.shape
    xm, xs = X.mean(0), X.std(0) + 1e-8
    ym, ys = Y.mean(), Y.std() + 1e-8
    Xs = torch.tensor((X - xm) / xs, dtype=torch.float32)
    Ys = torch.tensor((Y - ym) / ys, dtype=torch.float32)
    Tt = torch.tensor(T, dtype=torch.float32)

    # seeded, T-stratified validation split
    idx_t = np.flatnonzero(T == 1); idx_c = np.flatnonzero(T == 0)
    rng.shuffle(idx_t); rng.shuffle(idx_c)
    nv_t = int(round(val_frac * len(idx_t))); nv_c = int(round(val_frac * len(idx_c)))
    val_idx = np.concatenate([idx_t[:nv_t], idx_c[:nv_c]])
    tr_idx = np.concatenate([idx_t[nv_t:], idx_c[nv_c:]])
    tr_idx = torch.tensor(np.sort(tr_idx)); val_idx = torch.tensor(np.sort(val_idx))

    net = TarNet(p, propensity_head=dragon)
    opt = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=weight_decay)

    def loss_fn(idx):
        y0, y1, logit = net(Xs[idx])
        t = Tt[idx]
        yhat = t * y1 + (1 - t) * y0
        loss = F.mse_loss(yhat, Ys[idx])
        if dragon:
            loss = loss + alpha_prop * F.binary_cross_entropy_with_logits(logit, t)
        return loss

    best_val, best_state, best_epoch, bad = np.inf, None, 0, 0
    n_tr = len(tr_idx)
    g = torch.Generator().manual_seed(seed)
    for epoch in range(1, max_epochs + 1):
        net.train()
        perm = tr_idx[torch.randperm(n_tr, generator=g)]
        for s in range(0, n_tr, batch):
            b = perm[s:s + batch]
            opt.zero_grad()
            loss = loss_fn(b)
            loss.backward()
            opt.step()
        net.eval()
        with torch.no_grad():
            vl = float(loss_fn(val_idx))
        if vl < best_val - 1e-7:
            best_val, bad, best_epoch = vl, 0, epoch
            best_state = {k: v.detach().clone() for k, v in net.state_dict().items()}
        else:
            bad += 1
        if epoch >= min_epochs and bad >= patience:
            break
    if best_state is not None:
        net.load_state_dict(best_state)
    net.eval()
    with torch.no_grad():
        y0, y1, logit = net(Xs)
    mu0 = y0.numpy() * ys + ym
    mu1 = y1.numpy() * ys + ym
    e_hat = torch.sigmoid(logit).numpy() if dragon else None
    info = {"epochs": epoch, "best_epoch": best_epoch, "best_val": best_val}
    return mu0, mu1, e_hat, info
