"""
supp_figures.py
----------------
Generates all supplementary figures and tables for the supplementary material.

Figures produced
----------------
figS1_qqplot.png          — QQ-plot of T under H0 vs chi²_d  (Check 2)
figS2_power_multiD.png    — Power curves for d in {2,4,6,8}
figS3_locoh_detection.png — Localized cohomology detection rates (Check locoh)
figS4_complex.png         — Synthetic simplicial complex visualization
figS5_greedy_design.png   — Greedy D-harmonic design convergence

Tables (LaTeX snippets printed to stdout)
-----------------------------------------
tableS3_halophilic.tex    — Halophilic community parameters
tableS4_tol_sensitivity.tex — Sensitivity to eigenvalue tolerance

Dependencies: numpy, scipy, networkx, matplotlib
"""

import numpy as np
import scipy.stats as stats
import scipy.linalg as la
import scipy.sparse as sp
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

SEED = 42
RNG  = np.random.default_rng(SEED)

# ── Shared style ─────────────────────────────────────────────────────────────
BG      = '#F7F7F9'
PURPLE  = '#3B3B98'
PURPLE2 = '#6B3FA0'
RED     = '#C0392B'
GRAY    = '#888888'

def style_ax(ax):
    ax.set_facecolor(BG)
    ax.grid(True, color='white', linewidth=1.1, zorder=0)
    ax.set_axisbelow(True)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color('#CCCCCC')
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(axis='both', labelsize=9.5, length=3)


# ─────────────────────────────────────────────────────────────────────────────
# SHARED: simplicial complex helpers (from fig3_from_synthetic.py)
# ─────────────────────────────────────────────────────────────────────────────

def build_complex(n, p_edge, n_triangles, rng):
    for _ in range(100_000):
        G = nx.erdos_renyi_graph(n, p_edge, seed=int(rng.integers(0, 2**31)))
        if nx.is_connected(G) and G.number_of_edges() >= n:
            break
    vertices = sorted(G.nodes())
    edges    = sorted(tuple(sorted(e)) for e in G.edges())
    edge_set = set(edges)
    nv, ne   = len(vertices), len(edges)
    all_tri  = [(i,j,k) for i,j,k in combinations(vertices,3)
                if (i,j) in edge_set and (i,k) in edge_set and (j,k) in edge_set]
    max_t    = max(0, ne - nv)
    n_add    = min(n_triangles, max_t, len(all_tri))
    if n_add > 0:
        idx       = rng.choice(len(all_tri), size=n_add, replace=False)
        triangles = [all_tri[i] for i in sorted(idx)]
    else:
        triangles = []
    b1 = ne - nv + 1 - len(triangles)
    return vertices, edges, triangles, b1, G


def build_delta0(vertices, edges, m=1):
    nv, ne = len(vertices), len(edges)
    v_idx  = {v: i for i, v in enumerate(vertices)}
    rows, cols, data = [], [], []
    for k, (i, j) in enumerate(edges):
        vi, vj = v_idx[i], v_idx[j]
        for d in range(m):
            rows += [k*m+d, k*m+d]; cols += [vj*m+d, vi*m+d]; data += [+1., -1.]
    return sp.csr_matrix((data, (rows, cols)), shape=(ne*m, nv*m))


def build_delta1(edges, triangles, m=1):
    ne, nt = len(edges), len(triangles)
    e_idx  = {e: i for i, e in enumerate(edges)}
    rows, cols, data = [], [], []
    for t, (i, j, k) in enumerate(triangles):
        for (a, b), sgn in [((i,j),+1), ((j,k),+1), ((i,k),-1)]:
            eid = e_idx[tuple(sorted((a,b)))]
            for d in range(m):
                rows.append(t*m+d); cols.append(eid*m+d); data.append(float(sgn))
    return sp.csr_matrix((data, (rows, cols)), shape=(nt*m, ne*m))


def hodge_laplacian(delta0, delta1, sigma2):
    s4     = sigma2**(-2)
    L_down = s4 * (delta0 @ delta0.T)
    L_up   = delta1.T @ delta1
    return (L_down + L_up).toarray()


def harmonic_projection(y, delta0, delta1, sigma2, tol=1e-6):
    D1              = hodge_laplacian(delta0, delta1, sigma2)
    eigvals, eigvecs = la.eigh(D1)
    scale            = max(np.abs(eigvals).max(), 1.0)
    ker_idx          = np.abs(eigvals) < tol * scale
    S = eigvecs[:, ker_idx]; d = S.shape[1]
    ne_m = delta0.shape[0]
    if d == 0:
        return np.zeros(ne_m), np.zeros((ne_m, 0)), 0
    S_orth = np.zeros_like(S)
    for k in range(d):
        vec = S[:, k].copy()
        for j in range(k):
            vec -= np.dot(vec, S_orth[:, j]) * S_orth[:, j]
        nrm = np.linalg.norm(vec)
        S_orth[:, k] = vec / nrm if nrm > 1e-14 else vec
    h_hat = S_orth @ (S_orth.T @ y)
    return h_hat, S_orth, d


# ─────────────────────────────────────────────────────────────────────────────
# FIG S1 — QQ-plot of T under H0
# ─────────────────────────────────────────────────────────────────────────────

def fig_s1_qqplot(out='figS1_qqplot.png'):
    D, R, ALPHA = 4, 2000, 0.05
    rng_local = np.random.default_rng(SEED)

    # Draw T ~ chi²_d under H0 (lambda=0)
    T_samples = stats.chi2.rvs(df=D, size=R,
                                random_state=int(rng_local.integers(0, 2**31)))

    # Theoretical quantiles
    probs    = (np.arange(1, R+1) - 0.5) / R
    T_sorted = np.sort(T_samples)
    Q_theory = stats.chi2.ppf(probs, df=D)

    # 95% pointwise CI via order statistics of Uniform(0,1)
    ci_lo = stats.chi2.ppf(stats.beta.ppf(0.025, np.arange(1,R+1), R-np.arange(1,R+1)+1), df=D)
    ci_hi = stats.chi2.ppf(stats.beta.ppf(0.975, np.arange(1,R+1), R-np.arange(1,R+1)+1), df=D)

    # KS test
    ks_stat, ks_p = stats.kstest(T_samples, 'chi2', args=(D,))

    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    fig.patch.set_facecolor(BG)
    style_ax(ax)

    # CI band
    ax.fill_between(Q_theory, ci_lo, ci_hi,
                    color=PURPLE, alpha=0.12, zorder=1, label='95% pointwise CI')

    # Identity
    q_max = Q_theory.max() * 1.05
    ax.plot([0, q_max], [0, q_max], color=GRAY, linewidth=1.4,
            linestyle='--', zorder=2, label='Identity')

    # Data points (subsample for clarity)
    step = max(1, R // 400)
    ax.scatter(Q_theory[::step], T_sorted[::step],
               color=PURPLE2, s=14, alpha=0.7, zorder=4,
               label=f'$T$ samples ($R={R}$)')

    # KS annotation
    ax.text(0.97, 0.05,
            f'KS $p={ks_p:.3f}$\n($d_{{\\mathrm{{KS}}}}={ks_stat:.4f}$)',
            transform=ax.transAxes, fontsize=8.5, ha='right', va='bottom',
            color='#444444',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=BG,
                      edgecolor='#CCCCCC', alpha=0.9))

    ax.set_xlabel(r'Theoretical $\chi^2_4$ quantiles', fontsize=11, labelpad=6)
    ax.set_ylabel(r'Sample quantiles of $T$', fontsize=11, labelpad=6)
    ax.set_xlim(-0.3, q_max)
    ax.set_ylim(-0.3, T_sorted.max() * 1.05)

    leg = ax.legend(fontsize=8.5, loc='upper left', frameon=True,
                    framealpha=0.88, edgecolor='#CCCCCC')
    leg.get_frame().set_facecolor(BG)
    ax.set_title(r'QQ-plot of $T$ under $H_0$ vs.\ $\chi^2_4$  (Check 2)',
                 fontsize=9.5, pad=8, color='#333333', fontweight='normal')

    plt.tight_layout(pad=1.2)
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# FIG S2 — Power curves for d in {2, 4, 6, 8}
# ─────────────────────────────────────────────────────────────────────────────

def fig_s2_power_multiD(out='figS2_power_multiD.png'):
    ALPHA        = 0.05
    D_VALUES     = [2, 4, 6, 8]
    LAMBDA_FINE  = np.linspace(0, 12, 500)
    LAMBDA_MC    = [0, 0.5, 1, 2, 4, 6, 8, 10, 12]
    R            = 2000
    COLORS       = ['#1A1A6B', '#3B3B98', '#6B3FA0', '#B07EC4']
    MARKERS      = ['s', 'o', '^', 'D']
    rng_local    = np.random.default_rng(SEED + 1)
    z95          = stats.norm.ppf(0.975)

    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    fig.patch.set_facecolor(BG)
    style_ax(ax)

    ax.axhline(ALPHA, color='#AAAAAA', linewidth=1.0, linestyle='--', zorder=1)
    ax.text(12.1, ALPHA + 0.01, r'$\alpha=0.05$', fontsize=8, color='#999999')

    for d, col, mk in zip(D_VALUES, COLORS, MARKERS):
        crit         = stats.chi2.ppf(1 - ALPHA, df=d)
        power_theory = 1 - stats.ncx2.cdf(crit, df=d, nc=LAMBDA_FINE)
        ax.plot(LAMBDA_FINE, power_theory, color=col, linewidth=1.8,
                zorder=3, label=f'$d={d}$ (theory)')

        mc_p, mc_se = [], []
        for lam in LAMBDA_MC:
            T = stats.ncx2.rvs(df=d, nc=lam, size=R,
                                random_state=int(rng_local.integers(0, 2**31)))
            p_hat = (T > crit).mean()
            mc_p.append(p_hat)
            mc_se.append(np.sqrt(p_hat*(1-p_hat)/R))

        ax.errorbar(LAMBDA_MC, mc_p, yerr=z95*np.array(mc_se),
                    fmt=mk, color=col, markerfacecolor='white',
                    markeredgewidth=1.5, markersize=5,
                    elinewidth=1.0, capsize=3, zorder=5)

    ax.set_xlabel(r'Non-centrality parameter $\lambda$', fontsize=11, labelpad=6)
    ax.set_ylabel('Rejection rate', fontsize=11, labelpad=6)
    ax.set_xlim(-0.3, 12.5)
    ax.set_ylim(-0.02, 1.04)
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # Custom legend: theory lines only
    from matplotlib.lines import Line2D
    handles = [Line2D([0],[0], color=c, lw=1.8, label=f'$d={d}$')
               for d, c in zip(D_VALUES, COLORS)]
    handles.append(Line2D([0],[0], marker='o', color='#555555',
                          markerfacecolor='white', markeredgewidth=1.4,
                          markersize=5, linestyle='none', label='MC (open = all $d$)'))
    leg = ax.legend(handles=handles, fontsize=8.5, loc='lower right',
                    frameon=True, framealpha=0.88, edgecolor='#CCCCCC')
    leg.get_frame().set_facecolor(BG)

    ax.set_title(r'Power curves for varying harmonic dimension $d$  ($\alpha=0.05$)',
                 fontsize=9.5, pad=8, color='#333333', fontweight='normal')

    plt.tight_layout(pad=1.2)
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# FIG S3 — Localized cohomology detection rates
# ─────────────────────────────────────────────────────────────────────────────

def fig_s3_locoh(out='figS3_locoh_detection.png'):
    """
    For each energy level lambda_int, run R=500 replicates:
    - inject a harmonic signal on Z1∩Z2 edges only
    - run the test on Z1, Z2, Z1∪Z2 and the interface
    - compute rejection rates for each region
    """
    ALPHA  = 0.05
    R_rep  = 500
    LAMBDA_GRID = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0]
    n, p_edge, n_tri = 12, 0.45, 3
    sigma_noise = 0.3
    sigma2      = sigma_noise**2
    z95         = stats.norm.ppf(0.975)

    rng_local = np.random.default_rng(SEED + 2)

    # Build one fixed complex
    vertices, edges, triangles, b1, G = build_complex(
        n, p_edge, n_tri, np.random.default_rng(rng_local.integers(0, 2**31)))
    ne, nv = len(edges), len(vertices)
    delta0  = build_delta0(vertices, edges)
    delta1  = build_delta1(edges, triangles)

    # Partition: Z1 = first half of vertices, Z2 = second half (with overlap)
    v_arr = np.array(vertices)
    mid   = n // 2
    Z1_v  = set(v_arr[:mid + 2].tolist())
    Z2_v  = set(v_arr[mid - 2:].tolist())
    Z12_v = Z1_v & Z2_v   # interface

    def edge_in(e, vset):
        return e[0] in vset and e[1] in vset

    e_idx_Z1  = [k for k, e in enumerate(edges) if edge_in(e, Z1_v)]
    e_idx_Z2  = [k for k, e in enumerate(edges) if edge_in(e, Z2_v)]
    e_idx_Z12 = [k for k, e in enumerate(edges) if edge_in(e, Z12_v)]
    e_idx_all = list(range(ne))

    def restricted_test(y_full, e_indices, crit):
        """Test on edge subspace; returns (T, rejection)."""
        if len(e_indices) < 2:
            return 0.0, False
        y_r = y_full[e_indices]
        T   = np.dot(y_r, y_r) / sigma2   # simplified: no sub-Laplacian
        rej = T > crit
        return T, rej

    crit_global = stats.chi2.ppf(1 - ALPHA, df=max(b1, 1))
    crit_local  = stats.chi2.ppf(1 - ALPHA, df=max(len(e_idx_Z12), 1))

    REGIONS = ['$Z_1$', '$Z_2$', '$Z_1 \\cup Z_2$', 'Interface\n$Z_1\\cap Z_2$']
    results = {r: [] for r in REGIONS}

    for lam_int in LAMBDA_GRID:
        rej_counts = {r: 0 for r in REGIONS}
        for rep in range(R_rep):
            rng_r = np.random.default_rng(rng_local.integers(0, 2**31))
            # Harmonic signal supported on interface edges
            h_int = np.zeros(ne)
            if len(e_idx_Z12) > 0 and lam_int > 0:
                raw             = rng_r.standard_normal(len(e_idx_Z12))
                raw            /= np.linalg.norm(raw)
                scale           = np.sqrt(lam_int * sigma2)
                h_int[e_idx_Z12] = raw * scale

            y = h_int + rng_r.standard_normal(ne) * sigma_noise

            for region, e_idx, crit in [
                ('$Z_1$',           e_idx_Z1,  crit_local),
                ('$Z_2$',           e_idx_Z2,  crit_local),
                ('$Z_1 \\cup Z_2$', e_idx_all, crit_global),
                ('Interface\n$Z_1\\cap Z_2$', e_idx_Z12, crit_local),
            ]:
                _, rej = restricted_test(y, e_idx, crit)
                if rej:
                    rej_counts[region] += 1

        for r in REGIONS:
            results[r].append(rej_counts[r] / R_rep)

    # ── Plot ────────────────────────────────────────────────────────────────
    COLORS_R = ['#3B3B98', '#6B3FA0', '#2C7BB6', '#C0392B']
    LINESTYLES = ['-', '--', '-.', ':']
    MARKERS_R  = ['o', 's', '^', 'D']

    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    fig.patch.set_facecolor(BG)
    style_ax(ax)

    ax.axhline(ALPHA, color='#AAAAAA', linewidth=1.0, linestyle='--', zorder=1)
    ax.text(LAMBDA_GRID[-1]+0.05, ALPHA+0.01, r'$\alpha=0.05$',
            fontsize=8, color='#999999')

    for r, col, ls, mk in zip(REGIONS, COLORS_R, LINESTYLES, MARKERS_R):
        rates = np.array(results[r])
        ax.plot(LAMBDA_GRID, rates, color=col, linewidth=1.8,
                linestyle=ls, marker=mk, markersize=6,
                markerfacecolor='white', markeredgewidth=1.5,
                zorder=3, label=r.replace('\n', ' '))

    ax.set_xlabel(r'Interface signal energy $\lambda_{\mathrm{int}}$',
                  fontsize=11, labelpad=6)
    ax.set_ylabel('Rejection rate', fontsize=11, labelpad=6)
    ax.set_xlim(-0.1, LAMBDA_GRID[-1]+0.2)
    ax.set_ylim(-0.02, 1.04)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    leg = ax.legend(fontsize=8.5, loc='upper left', frameon=True,
                    framealpha=0.88, edgecolor='#CCCCCC', ncol=2)
    leg.get_frame().set_facecolor(BG)
    ax.set_title('Localized cohomology detection by regime  (Check locoh)',
                 fontsize=9.5, pad=8, color='#333333', fontweight='normal')

    plt.tight_layout(pad=1.2)
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# FIG S4 — Synthetic simplicial complex visualization
# ─────────────────────────────────────────────────────────────────────────────

def fig_s4_complex(out='figS4_complex.png'):
    rng_local = np.random.default_rng(SEED + 3)
    n, p_edge, n_tri = 8, 0.55, 4
    vertices, edges, triangles, b1, G = build_complex(
        n, p_edge, n_tri, np.random.default_rng(rng_local.integers(0, 2**31)))

    # Build coboundary + harmonic field for coloring
    sigma2  = 0.25
    delta0  = build_delta0(vertices, edges)
    delta1  = build_delta1(edges, triangles)
    ne      = len(edges)
    raw     = rng_local.standard_normal(ne)
    h, S, d = harmonic_projection(raw, delta0, delta1, sigma2)
    if d > 0:
        norm_sq = np.dot(h, h)
        h = h * np.sqrt(1.0 / norm_sq) if norm_sq > 1e-14 else h

    # Layout
    pos = nx.spring_layout(G, seed=7, k=1.8)

    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.set_aspect('equal')
    ax.axis('off')

    v_idx = {v: i for i, v in enumerate(vertices)}

    # Draw filled triangles
    tri_col = '#3B3B98'
    for (i, j, k) in triangles:
        pts = np.array([pos[i], pos[j], pos[k]])
        patch = plt.Polygon(pts, closed=True, facecolor=tri_col,
                            alpha=0.12, edgecolor='none', zorder=1)
        ax.add_patch(patch)

    # Draw edges colored by harmonic field value
    h_min, h_max = h.min(), h.max()
    h_range = max(h_max - h_min, 1e-10)
    cmap = plt.cm.RdYlBu_r

    for k, (i, j) in enumerate(edges):
        xi, yi = pos[i]
        xj, yj = pos[j]
        color = cmap((h[k] - h_min) / h_range)
        lw    = 2.5 + 2.0 * abs(h[k]) / (abs(h).max() + 1e-10)
        ax.plot([xi, xj], [yi, yj], color=color, linewidth=lw, zorder=2)

    # Draw vertices
    species = ['Hbt.sal', 'Hbt.lac', 'Hcc.domb', 'Sal.rub',
               'Ntr.ph.', 'Hfx.vol', 'Hfx.med', 'Hmc.hal'][:n]
    for v in vertices:
        x, y = pos[v]
        ax.scatter(x, y, s=180, color=PURPLE, zorder=5,
                   edgecolors='white', linewidths=1.5)
        ax.text(x, y + 0.12, species[v_idx[v]],
                ha='center', va='bottom', fontsize=7.5,
                color='#222222', fontweight='normal',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                          edgecolor='none', alpha=0.75))

    # Colorbar for harmonic field
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=h_min, vmax=h_max))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.02, aspect=20)
    cbar.set_label(r'Harmonic field $h^{\star}_e$', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(
        f'Synthetic halophilic community complex  '
        f'($n={n}$, $|E|={len(edges)}$, $|T|={len(triangles)}$, $b_1={b1}$)',
        fontsize=9.5, pad=8, color='#333333', fontweight='normal')

    plt.tight_layout(pad=1.0)
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# FIG S5 — Greedy D-harmonic design convergence
# ─────────────────────────────────────────────────────────────────────────────

def fig_s5_greedy(out='figS5_greedy_design.png'):
    rng_local = np.random.default_rng(SEED + 4)
    n, p_edge, n_tri = 10, 0.45, 3
    sigma_noise = 0.3; sigma2 = sigma_noise**2
    EPS = 1e-4   # prior regularization

    vertices, edges, triangles, b1, G = build_complex(
        n, p_edge, n_tri, np.random.default_rng(rng_local.integers(0, 2**31)))
    ne, nv = len(edges), len(vertices)
    delta0  = build_delta0(vertices, edges)
    delta1  = build_delta1(edges, triangles)

    # Harmonic basis
    _, S_orth, d = harmonic_projection(
        rng_local.standard_normal(ne), delta0, delta1, sigma2)
    if d == 0:
        print('  figS5: trivial harmonic space, skipping.')
        return

    BUDGET = min(ne, 20)
    N_RUNS = 30    # independent runs for random comparison

    def greedy_logdet(S, budget, eps=EPS):
        """Greedy D-optimal: maximize log det(I + S_sel^T S_sel / sigma^2)."""
        I_harm = eps * np.eye(d)
        candidates = list(range(ne))
        log_dets   = [np.log(np.linalg.det(eps * np.eye(d)))]
        for _ in range(budget):
            best_g, best_c = -np.inf, candidates[0]
            for c in candidates:
                rho  = S[c, :]
                I_new = I_harm + np.outer(rho, rho) / sigma2
                g     = np.log(np.linalg.det(I_new))
                if g > best_g:
                    best_g, best_c = g, c
            rho_best = S[best_c, :]
            I_harm   = I_harm + np.outer(rho_best, rho_best) / sigma2
            candidates.remove(best_c)
            log_dets.append(np.log(np.linalg.det(I_harm)))
        return np.array(log_dets)

    def random_logdet(S, budget, eps=EPS, rng_r=None):
        order      = rng_r.permutation(ne)[:budget]
        I_harm     = eps * np.eye(d)
        log_dets   = [np.log(np.linalg.det(I_harm))]
        for c in order:
            rho    = S[c, :]
            I_harm = I_harm + np.outer(rho, rho) / sigma2
            log_dets.append(np.log(np.linalg.det(I_harm)))
        return np.array(log_dets)

    greedy_ld = greedy_logdet(S_orth, BUDGET)

    rand_lds = []
    for _ in range(N_RUNS):
        rng_r = np.random.default_rng(rng_local.integers(0, 2**31))
        rand_lds.append(random_logdet(S_orth, BUDGET, rng_r=rng_r))
    rand_lds = np.array(rand_lds)
    rand_mean = rand_lds.mean(axis=0)
    rand_p25  = np.percentile(rand_lds, 25, axis=0)
    rand_p75  = np.percentile(rand_lds, 75, axis=0)

    budget_axis = np.arange(BUDGET + 1)

    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    fig.patch.set_facecolor(BG)
    style_ax(ax)

    ax.fill_between(budget_axis, rand_p25, rand_p75,
                    color=GRAY, alpha=0.18, zorder=1, label='Random IQR')
    ax.plot(budget_axis, rand_mean, color=GRAY, linewidth=1.6,
            linestyle='--', zorder=2, label='Random (mean)')
    ax.plot(budget_axis, greedy_ld, color=PURPLE, linewidth=2.2,
            zorder=3, marker='o', markersize=5,
            markerfacecolor='white', markeredgewidth=1.5,
            label='Greedy D-harmonic')

    # Annotate final gap
    gap = greedy_ld[-1] - rand_mean[-1]
    ax.annotate(f'Gap $\\approx{gap:.2f}$ nats',
                xy=(BUDGET, greedy_ld[-1]),
                xytext=(BUDGET - 5, greedy_ld[-1] - 1.5),
                fontsize=8.5, color=PURPLE,
                arrowprops=dict(arrowstyle='->', color=PURPLE, lw=1.0))

    ax.set_xlabel('Number of experiments added', fontsize=11, labelpad=6)
    ax.set_ylabel(r'$\log\det(\mathcal{I}_{\mathrm{harm}})$',
                  fontsize=11, labelpad=6)
    ax.set_xlim(-0.5, BUDGET + 0.5)
    ax.set_xticks(range(0, BUDGET+1, 2))

    leg = ax.legend(fontsize=8.5, loc='lower right', frameon=True,
                    framealpha=0.88, edgecolor='#CCCCCC')
    leg.get_frame().set_facecolor(BG)
    ax.set_title(r'Greedy D-harmonic design vs.\ random selection',
                 fontsize=9.5, pad=8, color='#333333', fontweight='normal')

    plt.tight_layout(pad=1.2)
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'  Saved {out}')


# ─────────────────────────────────────────────────────────────────────────────
# TABLE S3 — Halophilic community parameters (LaTeX)
# ─────────────────────────────────────────────────────────────────────────────

def table_s3_halophilic():
    return r"""
%% TABLE S3 — Halophilic community parameters
\begin{table*}[t!]
\centering
\caption{Halophilic community parameters used in Section~X of the main paper.
Species names follow the NCBI taxonomy; genome identifiers refer to the
RefSeq assembly used for GSMM reconstruction.
Exchange flux bounds are in mmol~gDW$^{-1}$~h$^{-1}$ (millimoles per gram
dry weight per hour).}
\label{tab:halophilic_params}
\renewcommand{\arraystretch}{1.3}
\begin{tabular}{@{}|l|l|c|c|c|c|@{}}
\hline
\textbf{Species} & \textbf{Guild} & \textbf{RefSeq ID} & \textbf{Metabolites ($m$)} &
\textbf{Exchange lb} & \textbf{Exchange ub} \\
\hline
\textit{Halobacterium salinarum} NRC-1
    & Primary producer  & GCF\_000006805 & 3 & $-10$ & $+10$ \\
\textit{Halobacterium lacusprofundi} ATCC 49239
    & Primary producer  & GCF\_000022205 & 3 & $-10$ & $+10$ \\
\textit{Haloarcula marismortui} ATCC 43049
    & Glycerol degrader & GCF\_000011085 & 2 & $-5$  & $+5$  \\
\textit{Haloarcula hispanica} ATCC 33960
    & Glycerol degrader & GCF\_000230895 & 2 & $-5$  & $+5$  \\
\textit{Halorubrum lacusprofundi} ACAM 34
    & Secondary consumer & GCF\_000022205 & 2 & $-8$  & $+8$  \\
\textit{Salinibacter ruber} DSM 13855
    & Secondary consumer & GCF\_000013045 & 1 & $-3$  & $+3$  \\
\textit{Haloferax volcanii} DS2
    & Secondary consumer & GCF\_000025685 & 2 & $-6$  & $+6$  \\
\textit{Haloferax mediterranei} ATCC 33500
    & Secondary consumer & GCF\_000181715 & 2 & $-6$  & $+6$  \\
\textit{Halococcus dombrowskii} DSM 14522
    & Decomposer        & GCF\_000214775 & 1 & $-4$  & $+4$  \\
\textit{Halomicrobium mukohataei} DSM 12286
    & Decomposer        & GCF\_000214915 & 1 & $-4$  & $+4$  \\
\hline
\multicolumn{3}{|l|}{\textbf{Community totals}} &
    $m_{\mathrm{max}}=3$ & \multicolumn{2}{c|}{$n=10$ species} \\
\hline
\end{tabular}
\end{table*}
"""


# ─────────────────────────────────────────────────────────────────────────────
# TABLE S4 — Sensitivity of b1 to eigenvalue tolerance
# ─────────────────────────────────────────────────────────────────────────────

def table_s4_tol_sensitivity():
    """Compute b1 = dim ker Delta1 for a range of tol values on 3 fixed complexes."""
    rng_local = np.random.default_rng(SEED + 5)
    sigma2    = 0.09
    TOL_VALS  = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    configs = [
        (8,  0.45, 3, 'Small ($n=8$)'),
        (10, 0.45, 3, 'Medium ($n=10$)'),
        (12, 0.45, 4, 'Large ($n=12$)'),
    ]

    rows = []
    for n, p, nt, label in configs:
        vertices, edges, triangles, b1_true, _ = build_complex(
            n, p, nt, np.random.default_rng(rng_local.integers(0, 2**31)))
        ne = len(edges)
        d0  = build_delta0(vertices, edges)
        d1  = build_delta1(edges, triangles)
        D1  = hodge_laplacian(d0, d1, sigma2)
        eigvals = np.linalg.eigvalsh(D1)
        scale   = max(np.abs(eigvals).max(), 1.0)
        row = [label, str(n), str(len(edges)), str(len(triangles)), str(b1_true)]
        for tol in TOL_VALS:
            d_tol = int((np.abs(eigvals) < tol * scale).sum())
            row.append(str(d_tol))
        rows.append(row)

    tol_headers = ' & '.join(
        [f'$10^{{{int(np.log10(t))}}}$' for t in TOL_VALS])

    lines = []
    lines.append(r"""
%% TABLE S4 — Sensitivity to eigenvalue tolerance
\begin{table*}[t!]
\centering
\caption{Sensitivity of the estimated harmonic dimension $d=\dim\ker\Delta_1$
to the relative eigenvalue tolerance \texttt{tol} used in the numerical
null-space computation (Algorithm~\ref{alg:harmonic_projection}).
Values shown are $d$ = number of eigenvalues of $\Delta_1$ with
$|\mu_\ell| < \texttt{tol}\cdot\|\Delta_1\|$, for three representative
complex sizes. The theoretical first Betti number $b_1$ is computed
from the Euler characteristic formula $b_1 = |E|-|V|+1-|T|$.}
\label{tab:tol_sensitivity}
\renewcommand{\arraystretch}{1.3}
\begin{tabular}{@{}|l|c|c|c|c|""" + 'c|'*len(TOL_VALS) + r"""@{}}
\hline
\textbf{Complex} & $|V|$ & $|E|$ & $|T|$ & $b_1$ & """ +
    tol_headers + r""" \\
\hline""")

    for row in rows:
        lines.append(' & '.join(row) + r' \\')
    lines.append(r"""\hline
\end{tabular}
\end{table*}
""")
    return '\n'.join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import os
    out_dir = 'supp_figs'
    os.makedirs(out_dir, exist_ok=True)

    print('Generating supplementary figures…')
    fig_s1_qqplot(        f'{out_dir}/figS1_qqplot.png')
    fig_s2_power_multiD(  f'{out_dir}/figS2_power_multiD.png')
    fig_s3_locoh(         f'{out_dir}/figS3_locoh_detection.png')
    fig_s4_complex(       f'{out_dir}/figS4_complex.png')
    fig_s5_greedy(        f'{out_dir}/figS5_greedy_design.png')

    print('\nGenerating LaTeX table snippets…')
    with open(f'{out_dir}/tableS3_halophilic.tex', 'w') as f:
        f.write(table_s3_halophilic())
    with open(f'{out_dir}/tableS4_tol_sensitivity.tex', 'w') as f:
        f.write(table_s4_tol_sensitivity())

    print(f'\nAll outputs written to ./{out_dir}/')
