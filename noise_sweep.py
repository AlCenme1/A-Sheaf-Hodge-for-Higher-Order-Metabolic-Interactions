"""Arbitrary-CV noise model versus model-derived Sigma (Section VI-C / Response to R2).
Replaces the FVA-derived standard deviations by sigma = max(CV*|x|, 1e-3) and
recomputes the exact test for the triad and for the ten-genus community.
Reads results/triad_base.json and results/comm_base.json."""
from paths import RESULTS
import json, numpy as np
from scipy.stats import chi2

CVS = (0.05, 0.1, 0.2, 0.5)
FLOOR = 1e-3

def cv_test(contrib, cv):
    T, d = 0.0, 0
    for c in contrib.values():
        v = np.array(list(c["values"].values()), float)
        if len(v) < 2:
            continue
        s = np.maximum(cv * np.abs(v), FLOOR); w = 1 / s**2
        mu = (w * v).sum() / w.sum()
        T += float((w * (v - mu) ** 2).sum()); d += len(v) - 1
    return T, d, float(chi2.sf(T, d))

def sweep():
    out = {}
    for name, f in (("triad", "triad_base.json"), ("community", "comm_base.json")):
        if not (RESULTS / f).exists():
            print(f"{f} not found: run the corresponding notebook first; skipping {name}")
            continue
        contrib = json.load(open(RESULTS / f))["contrib"]
        out[name] = [[cv, *cv_test(contrib, cv)] for cv in CVS]
    return out

if __name__ == "__main__":
    out = sweep()
    for k, rows in out.items():
        for cv, T, d, p in rows:
            print(f"{k:10s} CV={cv:<5} T={T:12.2f} d={d:5d} p={p:.3g}")
    json.dump(out, open(RESULTS / "noise_sweep.json", "w"), indent=1)
