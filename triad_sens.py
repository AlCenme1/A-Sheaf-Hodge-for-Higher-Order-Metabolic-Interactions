from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json, time

import triad_analysis as ta
from triad_lib import solve, org_exchange, pair_transfers

out = {"sensitivity": [], "transfer": {}}
grid = [dict(frac=1.0), dict(frac=0.95), dict(glyc_D=2.0), dict(glyc_D=10.0),
        dict(supp=0.02), dict(supp=0.1)]
for g in grid:
    ta.FRAC = g.get("frac", 0.99)
    t = time.time()
    r = ta.run(glyc_D=g.get("glyc_D", 5.0), supp_rate=g.get("supp", 0.05), verbose=False)
    top = sorted(r["contrib"].items(), key=lambda kv: -kv[1]["T"])[:3]
    rec = dict(setting=g, T=r["T"], d=r["d"], p=r["p"], radius=r["radius"],
               top=[("|".join(k), v["T"]) for k, v in top], growth=r["growth"])
    out["sensitivity"].append(rec)
    print(rec, f"{time.time()-t:.0f}s", flush=True)
ta.FRAC = 0.99

# transfer operator: slope of S->H DHA transfer w.r.t. glycerol uptake by S
def point(orgs, glyc_D=5.0, glyc_medium=0.0):
    com, fl, g = solve(orgs, glyc_D=glyc_D, glyc_medium=glyc_medium)
    bases = ["glyc_e", "dha_e"]
    f = org_exchange(fl, orgs, bases)
    med = {b: -fl.get(f"EX_{b}", 0.0) for b in bases}
    x = pair_transfers(f, med, orgs, bases)
    return x[("S", "H", "dha_e")], -f["S"]["glyc_e"]

for name, orgs, kw0, key in [("SH", ("S", "H"), 0.0, "glyc_medium"),
                             ("DSH", ("D", "S", "H"), 5.0, "glyc_D")]:
    slopes = []
    for h in (0.25, 0.5):
        lo = max(kw0 - h, 0.0); hi = kw0 + h
        xl, ul = point(orgs, **{key: lo}); xh, uh = point(orgs, **{key: hi})
        slopes.append((xh - xl) / (uh - ul) if abs(uh - ul) > 1e-9 else float("nan"))
    out["transfer"][name] = slopes
    print(name, slopes, flush=True)
json.dump(out, open(str(RESULTS/"triad_sens.json"), "w"), indent=1, default=str)
print("done", flush=True)
