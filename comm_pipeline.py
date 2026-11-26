"""Coherence analysis of the ten-genus community model (Section VII-D) and
reusable stages for the positive control (Section VII-B)."""
from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json, time, itertools, pickle, collections, numpy as np

from comm_lib import *
from comm_coords import transfers
from scipy.stats import chi2

EXCLUDE = {f"M_{x}_e" for x in INORG} | {"M_h2o_e", "M_h_e"}
FLOOR = 1e-3
FRAC_FVA = 0.99

def solve_contexts(C, contexts, extra_for=None, log=None):
    """extra_for: function(U) -> list of extra reactions or None."""
    res = {}; t0 = time.time()
    for n, U in enumerate(contexts):
        ctx = C.build(list(U), extra=extra_for(U) if extra_for else None)
        r = tradeoff(ctx)
        if r is None:
            continue
        vecs = C.secretion_matrix(ctx, [m for m in C.ext if m not in EXCLUDE])
        sec = {k: x for k, x in org_secretion(ctx, r["v"], vecs).items() if abs(x) > 1e-7}
        med = {m: -r["v"][ctx["ridx"][C.ex_rxn[m]]] for m in C.ext if m not in EXCLUDE}
        med = {m: x for m, x in med.items() if abs(x) > 1e-7}
        res[U] = dict(G=r["G"], mu=r["mu"], growth=r["growth"], sec=sec, medium=med, v=r["v"])
        if log and n % 20 == 0:
            print(f"  {n}/{len(contexts)} {U} {time.time()-t0:.0f}s", flush=True)
    return res

def coordinates(res, order):
    X = {U: transfers(res[U], U, order) for U in res}
    coords = collections.defaultdict(dict)
    for U, x in X.items():
        for k, v in x.items():
            coords[k][U] = v
    return coords

def needed_fva(res, coords):
    need = collections.defaultdict(set)
    for (i, j, m), d in coords.items():
        ctxs = [U for U in res if i in U and j in U]
        vals = [d.get(U, 0.0) for U in ctxs]
        if max(vals) - min(vals) > 1e-6:
            for U in ctxs:
                need[U].add((i, m)); need[U].add((j, m))
    return need

def run_fva(C, res, need, extra_for=None, log=None):
    W = {}; t0 = time.time()
    for n, (U, keys) in enumerate(need.items()):
        ctx = C.build(list(U), extra=extra_for(U) if extra_for else None)
        mets = sorted({m for _, m in keys})
        allvecs = C.secretion_matrix(ctx, mets)
        vecs = {k: allvecs[k] for k in keys if k in allvecs}
        fv = fva_tradeoff_pfba(ctx, res[U]["mu"], res[U]["v"], vecs, frac=FRAC_FVA)
        W[U] = {k: (hi - lo) for k, (lo, hi) in fv.items()}
        if log and n % 20 == 0:
            print(f"  FVA {n}/{len(need)} {U} {len(vecs)} obj {time.time()-t0:.0f}s", flush=True)
    return W

def coherence(res, coords, W):
    T = 0.0; d = 0; contrib = {}; radius = 0.0; radius_key = None
    for (i, j, m), dct in coords.items():
        ctxs = [U for U in res if i in U and j in U]
        if len(ctxs) < 2:
            continue
        v = np.array([dct.get(U, 0.0) for U in ctxs])
        s = np.array([max(0.5 * np.hypot(np.nan_to_num(W.get(U, {}).get((i, m), 0.0)),
                                         np.nan_to_num(W.get(U, {}).get((j, m), 0.0))), FLOOR) for U in ctxs])
        w = 1 / s**2; mu = (w * v).sum() / w.sum()
        t = float((w * (v - mu)**2).sum())
        T += t; d += len(ctxs) - 1
        spread = float(v.max() - v.min())
        if spread > radius:
            radius, radius_key = spread, (i, j, m)
        contrib[(i, j, m)] = dict(T=t, spread=spread,
                                  values={"+".join(U): round(float(x), 5) for U, x in zip(ctxs, v)},
                                  sigma={"+".join(U): round(float(x), 5) for U, x in zip(ctxs, s)})
    p = float(chi2.sf(T, d)) if d else float("nan")
    return dict(T=T, d=d, p=p, radius=radius, radius_key=radius_key, contrib=contrib)

def summarize(out, k=10):
    top = sorted(out["contrib"].items(), key=lambda kv: -kv[1]["T"])[:k]
    return [("|".join(key), round(v["T"], 2), round(v["spread"], 4)) for key, v in top]

if __name__ == "__main__":
    t0 = time.time()
    C = Community(); C.set_medium(minimal_medium())
    contexts = [c for k in (2, 3) for c in itertools.combinations(C.orgs, k)]
    print("stage 1: trade-off + pFBA", flush=True)
    res = solve_contexts(C, contexts, log=True)
    pickle.dump(res, open(str(RESULTS/"comm_base_res.pkl"), "wb"))
    coords = coordinates(res, C.orgs)
    need = needed_fva(res, coords)
    print("coords", len(coords), "M", len({k[2] for k in coords}), "FVA objectives", sum(len(x) for x in need.values()), f"{time.time()-t0:.0f}s", flush=True)
    W = run_fva(C, res, need, log=True)
    pickle.dump(W, open(str(RESULTS/"comm_base_W.pkl"), "wb"))
    out = coherence(res, coords, W)
    print(f"T={out['T']:.2f} d={out['d']} p={out['p']:.3g} radius={out['radius']:.4f} at {out['radius_key']}", flush=True)
    for row in summarize(out):
        print(row, flush=True)
    json.dump({"T": out["T"], "d": out["d"], "p": out["p"], "radius": out["radius"],
               "radius_key": "|".join(out["radius_key"]) if out["radius_key"] else None,
               "n_coords": len(coords), "M": sorted({k[2] for k in coords}),
               "contrib": {"|".join(k): v for k, v in out["contrib"].items()},
               "growth": {"+".join(U): {o: float(g) for o, g in r["growth"].items()} for U, r in res.items()},
               "runtime_s": time.time() - t0},
              open(str(RESULTS/"comm_base.json"), "w"), indent=1)
    print("done", f"{time.time()-t0:.0f}s", flush=True)
