"""Positive control (Section VII-B): in every context that contains both
Haloquadratum (i, lactate producer) and Halobacterium (k), Halobacterium
takes up eps units of lactate. The pair (Chromohalobacter, Haloquadratum)
lactate transfer is then changed only in the triad {Chromo, Halobact., Haloq.}."""
from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json, pickle, time, numpy as np

from comm_pipeline import *
from scipy.stats import chi2

I, K, J, MET = "Haloquadratum", "Halobacterium", "Chromohalobacter", "M_lac__L_e"
TARGET = tuple(sorted((J, I))) + (MET,)          # ('Chromohalobacter','Haloquadratum','M_lac__L_e')
EPS = [float(x) for x in sys.argv[1:]] or [1.0, 3.0, 10.0]
import os
OUTF = os.environ.get("PC_OUT", str(RESULTS/"comm_pc.json"))
WITH_BASE = os.environ.get("PC_BASE", "1") == "1"

C = Community(); C.set_medium(minimal_medium())
base = pickle.load(open(str(RESULTS/"comm_base_res.pkl"), "rb"))
Wb = pickle.load(open(str(RESULTS/"comm_base_W.pkl"), "rb"))
affected = [U for U in base if I in U and K in U]

def rank_of(out, key):
    order = sorted(out["contrib"], key=lambda k: -out["contrib"][k]["T"])
    return order.index(key) + 1 if key in order else None

def coord_p(out, key, res):
    c = out["contrib"].get(key)
    if c is None: return None
    n = sum(1 for U in res if key[0] in U and key[1] in U)
    return float(chi2.sf(c["T"], n - 1))

rows = []
for eps in ([0.0] if WITH_BASE else []) + EPS:
    t = time.time()
    extra = (lambda U: [(f"PC_LACSINK_{K}", {MET: -1.0}, eps, eps)] if (I in U and K in U) else None)
    res = dict(base); W = dict(Wb)
    if eps > 0:
        new = solve_contexts(C, affected, extra_for=extra)
        dropped = [U for U in affected if U not in new]
        res.update(new)
        for U in dropped: res.pop(U, None)
        coords = coordinates(res, C.orgs)
        need = needed_fva(res, coords)
        # FVA only where something changed: affected contexts, and new keys elsewhere
        need_new = {U: keys for U, keys in need.items()
                    if U in affected or not keys <= set(Wb.get(U, {}).keys())}
        need_new = {U: (keys if U in affected else keys - set(Wb.get(U, {}).keys())) for U, keys in need_new.items()}
        Wn = run_fva(C, res, need_new, extra_for=extra)
        for U, d in Wn.items():
            W[U] = {**(W.get(U, {}) if U not in affected else {}), **d}
    else:
        dropped = []
        coords = coordinates(res, C.orgs)
    out = coherence(res, coords, W)
    tc = out["contrib"].get(TARGET, {})
    row = dict(eps=eps, T=out["T"], d=out["d"], p=out["p"], radius=out["radius"],
               target_T=tc.get("T"), target_p=coord_p(out, TARGET, res), target_rank=rank_of(out, TARGET),
               target_values=tc.get("values"), target_sigma=tc.get("sigma"),
               top=summarize(out, 5), dropped=["+".join(U) for U in dropped], secs=time.time() - t)
    rows.append(row)
    print(json.dumps({k: row[k] for k in ("eps", "T", "d", "p", "target_T", "target_p", "target_rank", "dropped", "secs")}, default=str), flush=True)
    print("   top:", row["top"], flush=True)
json.dump(rows, open(OUTF, "w"), indent=1, default=str)
print("done", flush=True)
