"""
nerve_lib.py -- core library for the design-nerve (delta0) framework.
"""
import numpy as np
from itertools import combinations


def pairs_of(species_set):
    return sorted(combinations(sorted(species_set), 2))


class Context:
    __slots__ = ("species", "pairs", "pair_index")
    def __init__(self, species):
        self.species = tuple(sorted(species))
        self.pairs = pairs_of(self.species)
        self.pair_index = {p: k for k, p in enumerate(self.pairs)}
    def dim(self):
        return len(self.pairs)


def build_nerve(contexts):
    edges = []
    for u in range(len(contexts)):
        for v in range(u + 1, len(contexts)):
            shared = sorted(set(contexts[u].pairs) & set(contexts[v].pairs))
            if shared:
                edges.append((u, v, shared))
    return edges


def stalk_offsets(contexts):
    offs = [0]
    for c in contexts:
        offs.append(offs[-1] + c.dim())
    return offs


def build_delta0(contexts, edges):
    offs = stalk_offsets(contexts)
    n0 = offs[-1]
    e_offs = [0]
    for (u, v, shared) in edges:
        e_offs.append(e_offs[-1] + len(shared))
    n1 = e_offs[-1]
    D = np.zeros((n1, n0))
    for k, (u, v, shared) in enumerate(edges):
        row0 = e_offs[k]
        for r, p in enumerate(shared):
            D[row0 + r, offs[u] + contexts[u].pair_index[p]] += 1.0
            D[row0 + r, offs[v] + contexts[v].pair_index[p]] -= 1.0
    return D, offs, e_offs


def _pinv_quad(M, tol_rel=1e-9):
    w, V = np.linalg.eigh(M)
    tol = tol_rel * max(w.max(), 1e-30)
    keep = w > tol
    Minv = (V[:, keep] * (1.0 / w[keep])) @ V[:, keep].T
    return Minv, int(keep.sum())


def gls_test(delta0, Sigma_diag, y):
    Dy = delta0 @ y
    Gram = delta0 @ np.diag(Sigma_diag) @ delta0.T
    Ginv, d = _pinv_quad(Gram)
    return float(Dy @ Ginv @ Dy), d


def gls_test_indexed(delta0, Sigma, y, idx):
    """chi2 stat using only the edge-coordinates in idx."""
    A = delta0[idx, :]
    Ay = A @ y
    M = A @ Sigma @ A.T
    Minv, d = _pinv_quad(M)
    return float(Ay @ Minv @ Ay), d


def nested_regime_decomposition(delta0, Sigma, y, local_idx, if_idx):
    """Exact T_full = T_loc + T_(if|loc) via Schur complement."""
    A1 = delta0[local_idx, :]
    A2 = delta0[if_idx, :]
    M11 = A1 @ Sigma @ A1.T
    M12 = A1 @ Sigma @ A2.T
    M22 = A2 @ Sigma @ A2.T
    M11_inv, d_loc = _pinv_quad(M11)
    Ay1, Ay2 = A1 @ y, A2 @ y
    T_loc = float(Ay1 @ M11_inv @ Ay1)
    r = Ay2 - M12.T @ (M11_inv @ Ay1)
    Schur = M22 - M12.T @ M11_inv @ M12
    Schur_inv, d_cond = _pinv_quad(Schur)
    T_cond = float(r @ Schur_inv @ r)
    return T_loc, d_loc, T_cond, d_cond


def h0_projection(delta0, Sigma_diag, y):
    Sig = np.diag(Sigma_diag)
    Dy = delta0 @ y
    Gram = delta0 @ Sig @ delta0.T
    Ginv = np.linalg.pinv(Gram)
    return y - Sig @ delta0.T @ (Ginv @ Dy)


def orth_complement_vector(delta0, w_edge):
    return delta0.T @ w_edge
