from paths import FIGS
from plotstyle import NAVY, RUST, TEAL, GOLD, GRAY, LIGHT, INK, style_ax, savefig
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import json, pickle, itertools
from nerve_lib import Context, build_nerve, build_delta0, gls_test, orth_complement_vector, nested_regime_decomposition

# ============================================================ FIG 2 (main): power curve
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

rng = np.random.default_rng(1)
delta0, offs, e_offs = make_nerve(6, 8, 5, rng)
n0 = offs[-1]; sigma = 0.3
Sigma_diag = np.full(n0, sigma**2)
_, d_fixed = gls_test(delta0, Sigma_diag, np.zeros(n0))

lam_grid = [0, 0.5, 1, 2, 4, 8]; R = 800
emp = []
for lam in lam_grid:
    if lam > 0:
        w = rng.normal(size=delta0.shape[0])
        r = orth_complement_vector(delta0, w)
        T0, _ = gls_test(delta0, Sigma_diag, r)
        r = r * np.sqrt(lam / T0) if T0 > 1e-12 else r
    else:
        r = np.zeros(n0)
    crit = stats.chi2.ppf(0.95, d_fixed)
    rej = 0
    for _ in range(R):
        eps = rng.normal(scale=np.sqrt(Sigma_diag))
        T, _ = gls_test(delta0, Sigma_diag, r + eps)
        rej += T > crit
    emp.append(rej / R)
emp = np.array(emp)
lam_fine = np.linspace(0, max(lam_grid), 200)
theo = 1 - stats.ncx2.cdf(stats.chi2.ppf(0.95, d_fixed), d_fixed, lam_fine)

fig, ax = plt.subplots(figsize=(3.6, 2.9))
style_ax(ax)
ax.plot(lam_fine, theo, color=NAVY, lw=1.8, label=r'theoretical: $1-F_{\chi^2_d(\lambda)}(\chi^2_{d,0.95})$')
se = np.sqrt(emp * (1 - emp) / R)
ax.errorbar(lam_grid, emp, yerr=1.96*se, fmt='o', color=RUST, ecolor=RUST, ms=4.5,
            elinewidth=1, capsize=2.5, label=f'Monte Carlo ($R={R}$)')
ax.axhline(0.05, color=GRAY, ls=(0,(4,2)), lw=1, label=r'$\alpha=0.05$')
ax.set_xlabel(r'non-centrality $\lambda$')
ax.set_ylabel('rejection rate')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=2, frameon=False, handlelength=1.6)
savefig(fig, str(FIGS/'FIG2_power_curve.png'))
print("FIG2 done")

# ============================================================ FIG 3 (main): factorization residual
def factorization_residual(Y):
    s = np.linalg.svd(Y, compute_uv=False)
    return float(np.sum(s[1:]**2))

rng = np.random.default_rng(7)
p, q, n_rep = 6, 6, 80
R_targets = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
noise_sd = 0.04
null_vals = []
for _ in range(300):
    a1 = rng.normal(size=p); a2 = rng.normal(size=q)
    pure = np.outer(a1, a2); pure /= np.linalg.norm(pure)
    eps = rng.normal(scale=noise_sd, size=(p, q))
    null_vals.append(factorization_residual(pure + eps))
null_floor = np.mean(null_vals)

recovered_mean, recovered_ci = [], []
for Rt in R_targets:
    vals = []
    for _ in range(n_rep):
        a1 = rng.normal(size=p); a2 = rng.normal(size=q)
        pure = np.outer(a1, a2); pure /= np.linalg.norm(pure)
        u1, s1, v1t = np.linalg.svd(pure)
        u_top, v_top = u1[:, :1], v1t[:1, :]
        pert = rng.normal(size=(p, q))
        pert -= u_top @ (u_top.T @ pert)
        pert -= (pert @ v_top.T) @ v_top
        if np.linalg.norm(pert) > 1e-12:
            pert = pert / np.linalg.norm(pert) * np.sqrt(Rt)
        Y = pure + pert
        eps = rng.normal(scale=noise_sd, size=(p, q))
        vals.append(factorization_residual(Y + eps))
    vals = np.array(vals)
    recovered_mean.append(vals.mean())
    recovered_ci.append(1.96 * vals.std(ddof=1) / np.sqrt(n_rep))
recovered_mean = np.array(recovered_mean); recovered_ci = np.array(recovered_ci)

fig, ax = plt.subplots(figsize=(3.6, 2.9))
style_ax(ax)
ax.plot([0, 1.05], [null_floor, null_floor + 1.05], ls=(0,(4,2)), color=GRAY,
        lw=1, label=f'identity $+$ empirical null floor ({null_floor:.3f})')
ax.errorbar(R_targets, recovered_mean, yerr=recovered_ci, fmt='o', color=TEAL,
            ms=4.5, elinewidth=1, capsize=2.5, label=r'recovered $\widehat R$ (mean$\pm$95\% CI)')
ax.set_xlabel(r'injected residual $\|r\|_F^2$')
ax.set_ylabel(r'recovered $\widehat R$')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=2, frameon=False, handlelength=1.6)
savefig(fig, str(FIGS/'FIG3_factorization_residual.png'))
print("FIG3 done")

# ============================================================ FIG 4 (main): exact regime decomposition
rng = np.random.default_rng(11)
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
for (u, v, shared) in edges:
    ru, rv = ctx_regime(contexts[u]), ctx_regime(contexts[v])
    edge_regime.append('IF' if (ru != rv or ru == 'IF') else ru)
local_idx, if_idx = [], []
for k, (u, v, shared) in enumerate(edges):
    r0, r1 = e_offs[k], e_offs[k+1]
    (local_idx if edge_regime[k] in ('Z1','Z2') else if_idx).extend(range(r0, r1))

Sigma_diag = np.full(n0, 0.25**2); Sigma = np.diag(Sigma_diag)
w_if = np.zeros(e_offs[-1])
rng2 = np.random.default_rng(21)
w_if[if_idx] = rng2.normal(size=len(if_idx))
r = delta0.T @ w_if
T0, _ = gls_test(delta0, Sigma_diag, r)
r = r * np.sqrt(6.0 / T0)
eps = rng2.normal(scale=np.sqrt(Sigma_diag))
y = r + eps
T_loc, d_loc, T_cond, d_cond = nested_regime_decomposition(delta0, Sigma, y, local_idx, if_idx)
T_full, d_full = gls_test(delta0, Sigma_diag, y)
crit_loc = stats.chi2.ppf(0.95, d_loc); crit_cond = stats.chi2.ppf(0.95, d_cond); crit_full = stats.chi2.ppf(0.95, d_full)

fig, ax = plt.subplots(figsize=(3.6, 2.9))
style_ax(ax)
labels = ['$T_{\\mathrm{loc}}$', '$T_{\\mathrm{if|loc}}$', '$T_{\\mathrm{full}}$']
vals = [T_loc, T_cond, T_full]
crits = [crit_loc, crit_cond, crit_full]
colors = [NAVY, RUST, INK]
bars = ax.bar(labels, vals, color=colors, width=0.5, zorder=3)
for i, (v, c) in enumerate(zip(vals, crits)):
    ax.plot([i-0.28, i+0.28], [c, c], color=GRAY, ls=(0,(3,2)), lw=1.3, zorder=4)
    ax.text(i, v + 0.6, f'{v:.2f}', ha='center', fontsize=9.5)
ax.plot([], [], color=GRAY, ls=(0,(3,2)), lw=1.3, label=r'$\chi^2_{d,0.95}$ critical value')
ax.set_ylabel('test statistic')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=2, frameon=False, handlelength=1.6)
diff = T_full - (T_loc + T_cond)
ax.set_xlabel(f'$T_{{full}}-(T_{{loc}}+T_{{if|loc}})={diff:.1e}$', fontsize=8.5)
savefig(fig, str(FIGS/'FIG4_regime_decomposition.png'))
print("FIG4 done, diff=", diff)


# ============================================================ FIG 5 (main): real design comparison
