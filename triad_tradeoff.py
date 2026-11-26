"""Robustness: glycerol-DHA motif under the cooperative trade-off (alpha=0.5)
used for the community model, instead of maximal total growth."""
from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json; 
import triad_lib as tl, triad_analysis as ta
from cobra.flux_analysis import pfba as _pfba

def solve_tradeoff(orgs, alpha=0.5, **kw):
    com = tl.build_context(orgs, **kw)
    G = com.slim_optimize()
    gs = [com.reactions.get_by_id(f"Growth__{o}") for o in orgs if o != "D"]
    with com:
        for g in gs: g.lower_bound = alpha * G / len(gs)
        feasible = com.slim_optimize(error_value=None) is not None
    if feasible:
        for g in gs: g.lower_bound = 0.999 * alpha * G / len(gs)
    else:  # fall back: share the target among members that can grow
        for g in gs: g.lower_bound = 0.0
    sol = _pfba(com, fraction_of_optimum=alpha if not feasible else 0.0)
    fl = sol.fluxes.to_dict()
    return com, fl, {o: fl.get(f"Growth__{o}", 0.0) for o in orgs if o != "D"}

tl.solve = solve_tradeoff
ta.solve = solve_tradeoff
out = ta.run()
top = sorted(out["contrib"].items(), key=lambda kv: -kv[1]["T"])[:5]
print(f"T={out['T']:.2f} d={out['d']} p={out['p']:.3g} radius={out['radius']:.3f}")
for k, v in top: print(k, round(v["T"], 2), v["values"])
json.dump({"T": out["T"], "d": out["d"], "p": out["p"], "radius": out["radius"], "growth": out["growth"],
           "top": [("|".join(k), v["T"], v["values"]) for k, v in top]},
          open(str(RESULTS/"triad_tradeoff.json"), "w"), indent=1, default=float)
