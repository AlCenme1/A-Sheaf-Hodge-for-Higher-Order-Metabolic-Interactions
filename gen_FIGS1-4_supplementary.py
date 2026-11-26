from paths import FIGS
from plotstyle import NAVY, RUST, TEAL, GOLD, GRAY, LIGHT, INK, style_ax, savefig
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import itertools
from nerve_lib import Context, build_nerve, build_delta0, gls_test, orth_complement_vector, nested_regime_decomposition

def make_nerve(n_species, n_pair, n_triad, rng):
    species = [f"S{i}" for i in range(n_species)]
    pair_pool = list(itertools.combinations(range(n_species), 2))
    rng.shuffle(pair_pool)
    contexts = [Context([species[a], species[b]]) for a, b in pair_pool[:n_pair]]
    all_triples = list(itertools.combinations(range(n_species), 3))
    rng.shuffle(all_triples)
    for t in all_triples[:n_triad]:
        contexts.append(Context([species[i] for i in t]))
    edges = build_nerve(contexts)
    delta0, offs, e_offs = build_delta0(contexts, edges)
    return delta0, offs, e_offs

# ================================================================ FIG S1: QQ plot
rng = np.random.default_rng(2)
delta0, offs, e_offs = make_nerve(6, 8, 5, rng)
n0 = offs[-1]; Sigma_diag = np.full(n0, 0.3**2)
_, d = gls_test(delta0, Sigma_diag, np.zeros(n0))
R = 2000
Ts = np.empty(R)
for i in range(R):
    eps = rng.normal(scale=np.sqrt(Sigma_diag))
    Ts[i], _ = gls_test(delta0, Sigma_diag, eps)
Ts.sort()
theo_q = stats.chi2.ppf((np.arange(1, R+1)-0.5)/R, d)
ks_stat, ks_p = stats.kstest(Ts, 'chi2', args=(d,))

fig, ax = plt.subplots(figsize=(3.4, 3.1))
style_ax(ax, grid_axis='both')
lims = [0, max(Ts.max(), theo_q.max())*1.02]
ax.plot(lims, lims, ls=(0,(4,2)), color=GRAY, lw=1, label='identity')
ax.scatter(theo_q, Ts, s=8, color=TEAL, alpha=0.55, edgecolors='none', label=f'$T$ samples ($R={R}$)')
ax.set_xlabel(rf'theoretical $\chi^2_{{{d}}}$ quantiles')
ax.set_ylabel('sample quantiles of $T$')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), frameon=False, ncol=2, handlelength=1.6)
ax.text(0.96, 0.06, f'KS $p={ks_p:.3f}$', transform=ax.transAxes, ha='right', fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#DDDDDD', linewidth=0.7))
savefig(fig, str(FIGS/'FIGS1_qqplot.png'))
print("FIGS1 done, KS p=", ks_p)

# ================================================================ FIG S2: multi-d power
rng2 = np.random.default_rng(3)
configs = [(5,5,2), (6,8,5), (7,10,8), (8,14,10)]
lam_grid = [0,0.5,1,2,4,6,8,10,12]; R2 = 500
fig, ax = plt.subplots(figsize=(3.6, 2.95))
style_ax(ax)
colors = [NAVY, TEAL, GOLD, RUST]
for (n_sp, npair, ntri), col in zip(configs, colors):
    delta0, offs, e_offs = make_nerve(n_sp, npair, ntri, rng2)
    n0 = offs[-1]; Sigma_diag = np.full(n0, 0.3**2)
    _, d = gls_test(delta0, Sigma_diag, np.zeros(n0))
    emp = []
    for lam in lam_grid:
        if lam > 0:
            w = rng2.normal(size=delta0.shape[0])
            r = delta0.T @ w
            T0, _ = gls_test(delta0, Sigma_diag, r)
            r = r*np.sqrt(lam/T0) if T0 > 1e-12 else r
        else:
            r = np.zeros(n0)
        crit = stats.chi2.ppf(0.95, d); rej = 0
        for _ in range(R2):
            eps = rng2.normal(scale=np.sqrt(Sigma_diag))
            T, _ = gls_test(delta0, Sigma_diag, r+eps)
            rej += T > crit
        emp.append(rej/R2)
    emp = np.array(emp)
    lam_fine = np.linspace(0, max(lam_grid), 200)
    theo = 1 - stats.ncx2.cdf(stats.chi2.ppf(0.95, d), d, lam_fine)
    ax.plot(lam_fine, theo, color=col, lw=1.6, label=f'$d={d}$')
    se = np.sqrt(emp*(1-emp)/R2)
    ax.errorbar(lam_grid, emp, yerr=1.96*se, fmt='o', color=col, ms=3.3, elinewidth=0.8, capsize=1.8)
ax.axhline(0.05, color=GRAY, ls=(0,(4,2)), lw=1, label=r'$\alpha=0.05$')
ax.set_xlabel(r'non-centrality $\lambda$')
ax.set_ylabel('rejection rate')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), frameon=False, ncol=2, handlelength=1.6)
savefig(fig, str(FIGS/'FIGS2_power_multid.png'))
print("FIGS2 done")

# ================================================================ FIG S3: localized detection
rng = np.random.default_rng(31)
species8 = [f"S{i}" for i in range(8)]
Z1s, Z2s = species8[:4], species8[4:]
contexts = [Context([Z1s[0],Z1s[1]]), Context([Z1s[1],Z1s[2]]), Context([Z1s[0],Z1s[1],Z1s[2]]),
            Context([Z2s[0],Z2s[1]]), Context([Z2s[1],Z2s[2]]), Context([Z2s[0],Z2s[1],Z2s[2]]),
            Context([Z1s[1],Z2s[0]]), Context([Z1s[1],Z1s[2],Z2s[0]])]
edges = build_nerve(contexts)
delta0, offs, e_offs = build_delta0(contexts, edges)
n0 = offs[-1]
def ctx_regime(c):
    sp = set(c.species)
    if sp <= set(Z1s): return 'Z1'
    if sp <= set(Z2s): return 'Z2'
    return 'IF'
edge_regime = []
for (u,v,shared) in edges:
    ru, rv = ctx_regime(contexts[u]), ctx_regime(contexts[v])
    edge_regime.append('IF' if (ru != rv or ru == 'IF') else ru)
local_idx, if_idx = [], []
for k,(u,v,shared) in enumerate(edges):
    r0,r1 = e_offs[k], e_offs[k+1]
    (local_idx if edge_regime[k] in ('Z1','Z2') else if_idx).extend(range(r0,r1))
Sigma_diag = np.full(n0, 0.25**2); Sigma = np.diag(Sigma_diag)
lam_grid = [0,1,2,4,6,8,10,12]; R = 1000
# interface-only residual: r in ker(A_loc), so delta0 r vanishes on every local edge
import scipy.linalg as sla
A_loc = delta0[local_idx, :]
N = sla.null_space(A_loc)
r_dir = N @ rng.normal(size=N.shape[1])
assert np.abs(A_loc @ r_dir).max() < 1e-9 and np.abs(delta0 @ r_dir).max() > 1e-6
T_dir, _ = gls_test(delta0, Sigma_diag, r_dir)
rej_loc, rej_cond, rej_full = [], [], []
for lam in lam_grid:
    r = r_dir * np.sqrt(lam / T_dir)
    c_loc=c_cond=c_full=0
    for _ in range(R):
        eps = rng.normal(scale=np.sqrt(Sigma_diag))
        y = r+eps
        T_loc,d_loc,T_cond,d_cond = nested_regime_decomposition(delta0, Sigma, y, local_idx, if_idx)
        T_full,d_full = gls_test(delta0, Sigma_diag, y)
        c_loc += T_loc > stats.chi2.ppf(0.95, d_loc) if d_loc>0 else 0
        c_cond += T_cond > stats.chi2.ppf(0.95, d_cond) if d_cond>0 else 0
        c_full += T_full > stats.chi2.ppf(0.95, d_full)
    rej_loc.append(c_loc/R); rej_cond.append(c_cond/R); rej_full.append(c_full/R)

fig, ax = plt.subplots(figsize=(3.6, 2.95))
style_ax(ax)
ax.plot(lam_grid, rej_loc, 'o-', color=TEAL, ms=4, lw=1.4, label=r'$T_{\mathrm{loc}}$')
ax.plot(lam_grid, rej_full, '^-', color=NAVY, ms=4, lw=1.4, label=r'$T_{\mathrm{full}}$')
ax.plot(lam_grid, rej_cond, 'D-', color=RUST, ms=4, lw=1.4, label=r'$T_{\mathrm{if|loc}}$')
ax.axhline(0.05, color=GRAY, ls=(0,(4,2)), lw=1, label=r'$\alpha=0.05$')
ax.set_xlabel(r'injected interface-only signal $\lambda_{\mathrm{int}}$')
ax.set_ylabel('rejection rate')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), frameon=False, ncol=2, handlelength=1.6)
savefig(fig, str(FIGS/'FIGS3_locoh_detection.png'))
print("FIGS3 done")

# ================================================================ FIG S4: synthetic greedy vs random
rng3 = np.random.default_rng(51)
n_species_pool = 10
species_pool = [f"S{i}" for i in range(n_species_pool)]
base_pairs = list(itertools.combinations(range(n_species_pool), 2))
rng3.shuffle(base_pairs)
base_contexts = [Context([species_pool[a], species_pool[b]]) for a,b in base_pairs[:15]]
all_triples = list(itertools.combinations(range(n_species_pool), 3))
rng3.shuffle(all_triples)
candidate_contexts = [Context([species_pool[i] for i in t]) for t in all_triples[:30]]

def info_gain(current, candidate, sigma2=0.25**2):
    trial = current + [candidate]
    edges = build_nerve(trial)
    if not edges: return 0.0
    delta0, offs, e_offs = build_delta0(trial, edges)
    # regularized D-criterion (epsilon = 1): f = sum log(1 + eig(delta0^T Sigma^-1 delta0))
    w = np.linalg.eigvalsh(delta0 @ delta0.T) / sigma2
    w = np.clip(w, 0, None)
    return float(np.sum(np.log1p(w)))

def logdet_of(contexts, sigma2=0.25**2):
    edges = build_nerve(contexts)
    if not edges: return 0.0
    delta0, offs, e_offs = build_delta0(contexts, edges)
    # regularized D-criterion (epsilon = 1): f = sum log(1 + eig(delta0^T Sigma^-1 delta0))
    w = np.linalg.eigvalsh(delta0 @ delta0.T) / sigma2
    w = np.clip(w, 0, None)
    return float(np.sum(np.log1p(w)))

BUDGET = 12
current = list(base_contexts); remaining = list(candidate_contexts); trace = [logdet_of(current)]
for step in range(BUDGET):
    best_gain, best_c = -1e18, None
    for c in remaining:
        g = info_gain(current, c)
        if g > best_gain: best_gain, best_c = g, c
    current.append(best_c); remaining.remove(best_c); trace.append(best_gain)

n_rep = 15
random_traces = []
for _ in range(n_rep):
    cs = list(candidate_contexts); rng3.shuffle(cs)
    cur = list(base_contexts); tr = [logdet_of(cur)]
    for c in cs[:BUDGET]:
        cur.append(c); tr.append(logdet_of(cur))
    random_traces.append(tr)
random_traces = np.array(random_traces)
random_mean = random_traces.mean(axis=0)
random_q25 = np.percentile(random_traces, 25, axis=0)
random_q75 = np.percentile(random_traces, 75, axis=0)

fig, ax = plt.subplots(figsize=(3.6, 2.95))
style_ax(ax)
xs = np.arange(BUDGET+1)
ax.fill_between(xs, random_q25, random_q75, color=GRAY, alpha=0.22, label='random IQR')
ax.plot(xs, random_mean, ls=(0,(4,2)), color=GRAY, lw=1.2, label='random (mean)')
ax.plot(xs, trace, 'o-', color=NAVY, ms=4, lw=1.5, label='greedy D-optimal')
gap = trace[-1] - random_mean[-1]

ax.set_xlabel('experiments added')
ax.set_ylabel(r'$f=\log\det(I+\mathcal{I})$')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), frameon=False, ncol=2, handlelength=1.6)
savefig(fig, str(FIGS/'FIGS4_greedy_vs_random_synthetic.png'))
print(f"FIGS4 done, gap={gap:.3f}")
