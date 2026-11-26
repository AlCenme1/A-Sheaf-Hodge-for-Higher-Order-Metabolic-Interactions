from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json, time, itertools, numpy as np

from nerve_lib import Context, build_nerve, build_delta0
from comm_lib import Community
C = Community(); orgs = C.orgs
ext = set(C.ext)
org_mets = {o: {s for r in C.org_rxns[o] for s in (set(C.reactions[r]["reactants"]) | set(C.reactions[r]["products"])) if s in ext} for o in orgs}
n = len(orgs)
pairs = [(i, j) for i, j in itertools.combinations(range(n), 2) if org_mets[orgs[i]] & org_mets[orgs[j]]]
ps = set(pairs)
cand, overlap = [], []
for t in itertools.combinations(range(n), 3):
    i, j, k = t
    if (i, j) in ps and (j, k) in ps and (i, k) in ps:
        sh = org_mets[orgs[i]] & org_mets[orgs[j]] & org_mets[orgs[k]]
        if sh: cand.append(t); overlap.append(len(sh))
print("pairs", len(pairs), "candidates", len(cand), flush=True)
ctx = lambda t: Context([orgs[i] for i in t])
base = [ctx(p) for p in pairs]

def crit(contexts):
    edges = build_nerve(contexts)
    if not edges: return 0.0, -1e9, 0
    d0, offs, eo = build_delta0(contexts, edges)
    w = np.linalg.eigvalsh(d0 @ d0.T)
    pos = w[w > 1e-9 * max(w.max(), 1e-30)]
    return float(np.sum(np.log1p(np.clip(w, 0, None)))), float(np.sum(np.log(pos))), len(pos)

B = 24; t0 = time.time()
def greedy(key):
    cur = list(base); rem = list(range(len(cand))); ch = []
    for _ in range(B):
        best = max(rem, key=lambda i: crit(cur + [ctx(cand[i])])[key])
        cur.append(ctx(cand[best])); rem.remove(best); ch.append(best)
    return ch, crit(cur)
g_reg, c_reg = greedy(0)
g_pd, c_pd = greedy(1)
heur = sorted(range(len(cand)), key=lambda i: overlap[i])[:B]
c_h = crit(base + [ctx(cand[i]) for i in heur])
rng = np.random.default_rng(99)
rnd = np.array([crit(base + [ctx(cand[i]) for i in rng.choice(len(cand), B, replace=False)]) [:2] for _ in range(1000)])
out = dict(greedy_reg=c_reg, greedy_pdet=c_pd, heuristic=c_h,
           random_reg_mean=float(rnd[:, 0].mean()), random_reg_sd=float(rnd[:, 0].std()),
           random_pdet_mean=float(rnd[:, 1].mean()), random_pdet_sd=float(rnd[:, 1].std()),
           frac_random_ge_greedy_reg=float(np.mean(rnd[:, 0] >= c_reg[0])),
           frac_random_ge_heur_reg=float(np.mean(rnd[:, 0] >= c_h[0])),
           frac_random_ge_greedy_pdet=float(np.mean(rnd[:, 1] >= c_pd[1])),
           frac_random_ge_heur_pdet=float(np.mean(rnd[:, 1] >= c_h[1])),
           greedy_reg_triads=[[orgs[i] for i in cand[k]] for k in g_reg], secs=time.time() - t0)
print(json.dumps({k: v for k, v in out.items() if k != "greedy_reg_triads"}, indent=1), flush=True)
json.dump(out, open(str(RESULTS/"design_rev.json"), "w"), indent=1)
