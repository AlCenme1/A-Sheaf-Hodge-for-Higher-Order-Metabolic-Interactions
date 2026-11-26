"""Glycerol-DHA motif with flux-bounded FVA (pfba_factor=1.1), as for the community."""
from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json; 
import triad_analysis as ta
from triad_lib import flux_variability_analysis
def fva_widths(com, orgs, mets):
    rids = [f"LINK_{b}__{o}" for o in orgs if o != "D" for b in mets if f"LINK_{b}__{o}" in com.reactions]
    if not rids: return {}
    f = flux_variability_analysis(com, reaction_list=rids, fraction_of_optimum=ta.FRAC, pfba_factor=1.1)
    return {r: float(f.loc[r, "maximum"] - f.loc[r, "minimum"]) for r in rids}
ta.fva_widths = fva_widths
out = {}
for frac in (0.99, 0.95):
    ta.FRAC = frac
    r = ta.run(verbose=False)
    top = sorted(r["contrib"].items(), key=lambda kv: -kv[1]["T"])[:4]
    out[frac] = dict(T=r["T"], d=r["d"], p=r["p"], radius=r["radius"], top=[("|".join(k), v["T"], v["values"], v["sigma"]) for k, v in top])
    print(frac, f"T={r['T']:.2f} d={r['d']} p={r['p']:.3g} r={r['radius']:.3f}", flush=True)
    for t in out[frac]["top"]: print("   ", t, flush=True)
json.dump(out, open(str(RESULTS/"triad_pfbafactor.json"), "w"), indent=1, default=float)
