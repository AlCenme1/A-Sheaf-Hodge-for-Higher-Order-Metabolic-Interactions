from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json, time, numpy as np

from triad_lib import *
from scipy.stats import chi2

ORDER = ("D", "S", "H")
CONTEXTS = [("D", "S"), ("D", "H"), ("S", "H"), ("D", "S", "H")]
FRAC = 0.99
FLOOR = 1e-3

def context_data(orgs, **kw):
    com, fl, g = solve(orgs, **kw)
    bases = sorted({r.id[len("LINK_"):].rsplit("__", 1)[0] for r in com.reactions if r.id.startswith("LINK_")})
    if "D" in orgs and "glyc_e" not in bases:
        bases.append("glyc_e")
    f = org_exchange(fl, orgs, bases)
    med = {b: -fl.get(f"EX_{b}", 0.0) for b in bases}
    x = pair_transfers(f, med, orgs, bases, ORDER)
    return com, fl, g, f, med, x, bases

def fva_widths(com, orgs, mets):
    rids = [f"LINK_{b}__{o}" for o in orgs if o != "D" for b in mets if f"LINK_{b}__{o}" in com.reactions]
    if not rids:
        return {}
    fva = flux_variability_analysis(com, reaction_list=rids, fraction_of_optimum=FRAC, loopless=False)
    return {r: float(fva.loc[r, "maximum"] - fva.loc[r, "minimum"]) for r in rids}

def run(glyc_D=5.0, supp_rate=0.05, do_fva=True, verbose=True):
    res = {}
    for orgs in CONTEXTS:
        t = time.time()
        com, fl, g, f, med, x, bases = context_data(orgs, glyc_D=glyc_D, supp_rate=supp_rate)
        res[orgs] = dict(com=com, g=g, f=f, x=x)
        if verbose:
            print(orgs, "growth", {k: round(v, 4) for k, v in g.items()}, f"{time.time()-t:.1f}s", flush=True)
    # metabolite set: any nonzero transfer in any context
    M = sorted({k[2] for o in res for k, v in res[o]["x"].items() if abs(v) > 1e-6})
    pairs = sorted({(k[0], k[1]) for o in res for k in res[o]["x"]}, key=lambda p: (ORDER.index(p[0]), ORDER.index(p[1])))
    # widths
    for orgs in CONTEXTS:
        res[orgs]["w"] = fva_widths(res[orgs]["com"], orgs, M) if do_fva else {}
    # assemble coordinates: (pair, m) -> list of (context, value, sigma)
    coords = {}
    for orgs in CONTEXTS:
        for (i, j, b), v in res[orgs]["x"].items():
            if b not in M:
                continue
            wi = res[orgs]["w"].get(f"LINK_{b}__{i}", 0.0)
            wj = res[orgs]["w"].get(f"LINK_{b}__{j}", 0.0)
            sig = max(0.5 * np.hypot(wi, wj), FLOOR)
            coords.setdefault((i, j, b), []).append(("".join(orgs), v, sig))
    # exact test: for diagonal Sigma, ker(delta0) = one common value per coordinate,
    # T = sum over coordinates of weighted between-context dispersion
    T = 0.0; d = 0; contrib = {}; radius = 0.0
    for key, obs in coords.items():
        if len(obs) < 2:
            continue
        v = np.array([o[1] for o in obs]); s = np.array([o[2] for o in obs])
        w = 1 / s**2; mu = (w * v).sum() / w.sum()
        t = float((w * (v - mu) ** 2).sum())
        T += t; d += len(obs) - 1
        contrib[key] = dict(T=t, values={o[0]: round(o[1], 4) for o in obs},
                            sigma={o[0]: round(o[2], 4) for o in obs})
        radius = max(radius, float(v.max() - v.min()))
    p = float(chi2.sf(T, d)) if d > 0 else float("nan")
    return dict(T=T, d=d, p=p, radius=radius, M=M, contrib=contrib,
                growth={"".join(o): res[o]["g"] for o in CONTEXTS})

if __name__ == "__main__":
    out = run()
    top = sorted(out["contrib"].items(), key=lambda kv: -kv[1]["T"])[:10]
    print(f"T={out['T']:.2f} d={out['d']} p={out['p']:.3g} radius={out['radius']:.4f} |M|={len(out['M'])}")
    for k, v in top:
        print(k, round(v["T"], 2), v["values"], v["sigma"])
    json.dump({"T": out["T"], "d": out["d"], "p": out["p"], "radius": out["radius"], "M": out["M"],
               "growth": out["growth"],
               "contrib": {"|".join(k): v for k, v in out["contrib"].items()}},
              open(str(RESULTS/"triad_base.json"), "w"), indent=1)
