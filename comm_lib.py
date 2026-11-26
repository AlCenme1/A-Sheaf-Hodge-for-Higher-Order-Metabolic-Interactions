"""Community-context LP engine (HiGHS) for the ten-genus halophilic model.
Contexts are subsets of organisms; they share the 'ext' pool and the 544
R_EX_ medium exchanges. Net secretion of m_e by organism o is the sum over
o's reactions of S[m_e, r] v_r, which is linear in v."""
from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import pickle, numpy as np, scipy.sparse as sp, highspy
from pathlib import Path

INORG = ["h2o","h","na1","k","cl","mg2","ca2","so4","pi","nh4","fe2","fe3","mn2",
         "zn2","cu2","cobalt2","mobd","ni2","o2","co2"]
TRACE = ["fe2","fe3","mn2","zn2","cu2","cobalt2","mobd","ni2"]
INF = 1e30
PFBA_FRAC = 0.999   # pFBA within 0.1% of optimal growth

def minimal_medium(glyc=10.0, trace=0.1):
    """Minimal hypersaline medium: macro-ions and gases unlimited, CO2 only
    secreted, trace metals capped, glycerol as the sole carbon source."""
    med = {x: 1000.0 for x in INORG}
    med["co2"] = 0.0
    for t in TRACE:
        med[t] = trace
    med["glyc"] = glyc
    return med

class Community:
    def __init__(self, pkl=None, sbml=None):
        pkl = Path(pkl or COMM_PKL)
        if not pkl.exists():                      # parse once and cache
            from parse_sbml import parse_sbml
            pickle.dump(parse_sbml(str(sbml or COMMUNITY)), open(pkl, "wb"))
        self.species, self.reactions, self.params = pickle.load(open(pkl, "rb"))
        sp_, rx = self.species, self.reactions
        self.orgs = sorted({r[len("Growth_"):] for r in rx if r.startswith("Growth_")})
        # The SBML marks the shared pool metabolites as boundary species; they are
        # balanced here (only the community biomass sink stays a boundary).
        self.ext = sorted(s for s, v in sp_.items() if v.get("compartment") == "ext" and s != "community_biomass")
        self.balanced = set(self.ext)
        self.ex_rxn = {}                     # ext met -> R_EX_ id
        for r in rx:
            if r.startswith("R_EX_"):
                mets = set(rx[r]["reactants"]) | set(rx[r]["products"])
                for m in mets:
                    if m in self.ext:
                        self.ex_rxn[m] = r
        self.org_rxns = {o: sorted(r for r in rx if r.endswith("_" + o) and not r.startswith("R_EX_"))
                         for o in self.orgs}
        self.medium = {}                     # R_EX_ id -> max uptake (>=0)

    def bound(self, ref, default):
        v = self.params.get(ref) if ref is not None else None
        return default if v is None else v

    def set_medium(self, uptake):
        """uptake: dict ext-met-base (e.g. 'glyc') -> max uptake rate."""
        self.medium = {}
        for base, rate in uptake.items():
            m = f"M_{base}_e"
            if m in self.ex_rxn:
                self.medium[self.ex_rxn[m]] = rate

    def build(self, orgs, extra=None):
        """Assemble the LP for a context. extra: list of (rid, stoich dict, lb, ub)
        added reactions (used for the positive-control mechanism)."""
        rx = self.reactions
        rids = [r for o in orgs for r in self.org_rxns[o]] + sorted(self.ex_rxn.values())
        rids = list(dict.fromkeys(rids))
        stoich = {r: (rx[r]["reactants"], rx[r]["products"]) for r in rids}
        lb, ub = [], []
        for r in rids:
            if r.startswith("R_EX_"):
                lb.append(-self.medium.get(r, 0.0)); ub.append(1000.0)
            else:
                lb.append(self.bound(rx[r]["lb_ref"], -1000.0)); ub.append(self.bound(rx[r]["ub_ref"], 1000.0))
        if extra:
            for rid, st, l, u in extra:
                rids.append(rid); stoich[rid] = ({k: -v for k, v in st.items() if v < 0}, {k: v for k, v in st.items() if v > 0})
                lb.append(l); ub.append(u)
        mets = sorted({s for r in rids for s in (set(stoich[r][0]) | set(stoich[r][1]))
                       if (not self.species.get(s, {}).get("boundary", False)) or s in self.balanced})
        mi = {m: k for k, m in enumerate(mets)}
        rows, cols, data = [], [], []
        for k, r in enumerate(rids):
            for s, c in stoich[r][0].items():
                if s in mi: rows.append(mi[s]); cols.append(k); data.append(-c)
            for s, c in stoich[r][1].items():
                if s in mi: rows.append(mi[s]); cols.append(k); data.append(c)
        S = sp.csc_matrix((data, (rows, cols)), shape=(len(mets), len(rids)))
        ctx = dict(orgs=tuple(orgs), rids=rids, ridx={r: k for k, r in enumerate(rids)}, mets=mets, mi=mi,
                   S=S, lb=np.array(lb), ub=np.array(ub))
        ctx["growth"] = [ctx["ridx"][f"Growth_{o}"] for o in orgs]
        return ctx

    def secretion_vector(self, ctx, o, m):
        """Row vector c with c.v = net secretion of ext metabolite m by organism o."""
        S, mi, ridx = ctx["S"], ctx["mi"], ctx["ridx"]
        c = np.zeros(len(ctx["rids"]))
        if m not in mi:
            return c
        row = S.getrow(mi[m]) if hasattr(S, "getrow") else S[mi[m], :]
        row = sp.csr_matrix(S)[mi[m], :]
        for k, v in zip(row.indices, row.data):
            r = ctx["rids"][k]
            if not r.startswith("R_EX_") and r.endswith("_" + o):
                c[k] = v
        return c

    def secretion_matrix(self, ctx, mets):
        """Dict (o, m) -> sparse vector; only nonzero entries kept."""
        Sr = sp.csr_matrix(ctx["S"]); out = {}
        for m in mets:
            if m not in ctx["mi"]:
                continue
            row = Sr[ctx["mi"][m], :]
            for o in ctx["orgs"]:
                ks = [(k, v) for k, v in zip(row.indices, row.data)
                      if ctx["rids"][k].endswith("_" + o) and not ctx["rids"][k].startswith("R_EX_")]
                if ks:
                    out[(o, m)] = ks
        return out

def _highs(ctx, cost, extra_row=None):
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    lp = highspy.HighsLp()
    n = len(ctx["rids"]); A = ctx["S"]
    lp.num_col_ = n; lp.num_row_ = A.shape[0]
    lp.col_cost_ = np.asarray(cost, float); lp.col_lower_ = ctx["lb"]; lp.col_upper_ = ctx["ub"]
    lp.row_lower_ = np.zeros(A.shape[0]); lp.row_upper_ = np.zeros(A.shape[0])
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    h.passModel(lp)
    return h

def fba(ctx):
    """Maximize total growth; then pFBA (min sum |v|) at that optimum."""
    n = len(ctx["rids"]); c = np.zeros(n); c[ctx["growth"]] = -1.0
    h = _highs(ctx, c); h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        return None
    opt = -h.getInfo().objective_function_value
    # pFBA: split v = p - q, p,q >= 0, with bounds from lb/ub
    A = ctx["S"]; lb, ub = ctx["lb"], ctx["ub"]
    plo = np.maximum(lb, 0); pup = np.maximum(ub, 0)
    qlo = np.maximum(-ub, 0); qup = np.maximum(-lb, 0)
    g = np.zeros(n); g[ctx["growth"]] = 1.0
    A2 = sp.vstack([sp.hstack([A, -A]), sp.csr_matrix(np.concatenate([g, -g]))]).tocsc()
    h2 = highspy.Highs(); h2.setOptionValue("output_flag", False)
    lp = highspy.HighsLp(); lp.num_col_ = 2 * n; lp.num_row_ = A2.shape[0]
    lp.col_cost_ = np.ones(2 * n); lp.col_lower_ = np.concatenate([plo, qlo]); lp.col_upper_ = np.concatenate([pup, qup])
    lp.row_lower_ = np.concatenate([np.zeros(A.shape[0]), [PFBA_FRAC * opt]])
    lp.row_upper_ = np.concatenate([np.zeros(A.shape[0]), [INF]])
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A2.indptr; lp.a_matrix_.index_ = A2.indices; lp.a_matrix_.value_ = A2.data
    h2.passModel(lp); h2.run()
    if h2.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        v = np.array(h.getSolution().col_value)
    else:
        x = np.array(h2.getSolution().col_value); v = x[:n] - x[n:]
    return dict(opt=opt, v=v, growth={o: v[k] for o, k in zip(ctx["orgs"], ctx["growth"])})

def fva_secretion(ctx, opt, vecs, frac=0.99):
    """vecs: dict key -> list of (col, coef). Returns key -> (min, max) of c.v
    subject to total growth >= frac*opt."""
    n = len(ctx["rids"])
    h = _highs(ctx, np.zeros(n))
    g_idx = np.array(ctx["growth"], dtype=np.int32)
    h.addRow(frac * opt, INF, len(g_idx), g_idx, np.ones(len(g_idx)))
    out = {}; prev = np.array([], dtype=np.int32)
    for key, ks in vecs.items():
        idx = np.array([k for k, _ in ks], dtype=np.int32); val = np.array([v for _, v in ks], float)
        if len(prev):
            h.changeColsCost(len(prev), prev, np.zeros(len(prev)))
        res = []
        for sgn in (1.0, -1.0):
            h.changeColsCost(len(idx), idx, sgn * val)
            h.run()
            ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
            res.append(sgn * h.getInfo().objective_function_value if ok else np.nan)
        out[key] = (res[0], res[1])
        prev = idx
    return out

def org_secretion(ctx, v, vecs):
    return {key: float(sum(c * v[k] for k, c in ks)) for key, ks in vecs.items()}

def _lp_from(ctx, cost, extra_rows=()):
    """HighsLp for the context with optional extra rows (lo, hi, idx, val)."""
    n = len(ctx["rids"]); A = sp.csr_matrix(ctx["S"])
    rows = [A]; rlo = [np.zeros(A.shape[0])]; rhi = [np.zeros(A.shape[0])]
    for lo, hi, idx, val in extra_rows:
        r = np.zeros(n); r[idx] = val
        rows.append(sp.csr_matrix(r)); rlo.append([lo]); rhi.append([hi])
    A2 = sp.vstack(rows).tocsc()
    lp = highspy.HighsLp(); lp.num_col_ = n; lp.num_row_ = A2.shape[0]
    lp.col_cost_ = np.asarray(cost, float); lp.col_lower_ = ctx["lb"]; lp.col_upper_ = ctx["ub"]
    lp.row_lower_ = np.concatenate(rlo); lp.row_upper_ = np.concatenate(rhi)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A2.indptr; lp.a_matrix_.index_ = A2.indices; lp.a_matrix_.value_ = A2.data
    return lp

def tradeoff(ctx, alpha=0.5, tol=1e-7, max_iter=150):
    """MICOM cooperative trade-off (Diener et al. 2020): G* = max sum mu;
    mu* = argmin ||mu||^2 s.t. sum mu >= alpha G*, mu attainable. The QP is
    solved by Frank-Wolfe over the attainable growth vectors, with a warm-
    started LP oracle (the objective involves only the growth columns).
    Then pFBA with mu_o >= PFBA_FRAC mu*_o."""
    n = len(ctx["rids"]); gi = np.array(ctx["growth"], dtype=np.int32)
    c = np.zeros(n); c[gi] = -1.0
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.passModel(_lp_from(ctx, c)); h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        return None
    G = -h.getInfo().objective_function_value
    mu = np.array(h.getSolution().col_value)[gi]
    k_ = len(gi); eq = alpha * G / k_
    # exact shortcut: if every member can reach alpha G*/k at once, the equal
    # split is the minimum-norm point of {sum mu >= alpha G*}
    he = highspy.Highs(); he.setOptionValue("output_flag", False)
    he.passModel(_lp_from(ctx, np.zeros(n), [(eq, INF, np.array([g]), np.array([1.0])) for g in gi])); he.run()
    it = -1; gap = 0.0
    if he.getModelStatus() == highspy.HighsModelStatus.kOptimal:
        mu = np.full(k_, eq); max_iter = 0
    ho = highspy.Highs(); ho.setOptionValue("output_flag", False)
    ho.passModel(_lp_from(ctx, np.zeros(n), [(alpha * G, INF, gi, np.ones(len(gi)))]))
    for it in range(max_iter):
        grad = 2 * mu
        ho.changeColsCost(len(gi), gi, grad); ho.run()
        s_ = np.array(ho.getSolution().col_value)[gi]
        d = s_ - mu
        gap = -(grad @ d)
        if gap < tol * max(1.0, mu @ mu) * 1e3:
            break
        gam = min(1.0, max(0.0, -(mu @ d) / (d @ d)))
        mu = mu + gam * d
    A = ctx["S"]; lb, ub = ctx["lb"], ctx["ub"]
    rows = [sp.hstack([A, -A])]; rlo = [np.zeros(A.shape[0])]; rhi = [np.zeros(A.shape[0])]
    for k, m in zip(gi, mu):
        r = np.zeros(2 * n); r[k] = 1.0; r[n + k] = -1.0
        rows.append(sp.csr_matrix(r)); rlo.append([PFBA_FRAC * max(m, 0.0)]); rhi.append([INF])
    A2 = sp.vstack(rows).tocsc()
    lp = highspy.HighsLp(); lp.num_col_ = 2 * n; lp.num_row_ = A2.shape[0]
    lp.col_cost_ = np.ones(2 * n)
    lp.col_lower_ = np.concatenate([np.maximum(lb, 0), np.maximum(-ub, 0)])
    lp.col_upper_ = np.concatenate([np.maximum(ub, 0), np.maximum(-lb, 0)])
    lp.row_lower_ = np.concatenate(rlo); lp.row_upper_ = np.concatenate(rhi)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A2.indptr; lp.a_matrix_.index_ = A2.indices; lp.a_matrix_.value_ = A2.data
    hp = highspy.Highs(); hp.setOptionValue("output_flag", False); hp.passModel(lp); hp.run()
    x = np.array(hp.getSolution().col_value); v = x[:n] - x[n:]
    return dict(G=G, mu=dict(zip(ctx["orgs"], mu)), v=v, fw_iter=it, fw_gap=gap,
                growth={o: v[k] for o, k in zip(ctx["orgs"], gi)})

def fva_tradeoff(ctx, mu, vecs, frac=0.99):
    """FVA of secretion functionals with every organism's growth >= frac*mu_o."""
    n = len(ctx["rids"]); gi = ctx["growth"]
    rows = [(frac * mu[o], INF, np.array([k]), np.array([1.0])) for o, k in zip(ctx["orgs"], gi)]
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.setOptionValue("simplex_strategy", 4)   # primal simplex: only costs change between objectives
    h.passModel(_lp_from(ctx, np.zeros(n), rows))
    out = {}; prev = np.array([], dtype=np.int32)
    for key, ks in vecs.items():
        idx = np.array([k for k, _ in ks], dtype=np.int32); val = np.array([v for _, v in ks], float)
        if len(prev):
            h.changeColsCost(len(prev), prev, np.zeros(len(prev)))
        res = []
        for sgn in (1.0, -1.0):
            h.changeColsCost(len(idx), idx, sgn * val); h.run()
            ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
            res.append(sgn * h.getInfo().objective_function_value if ok else np.nan)
        out[key] = (res[0], res[1]); prev = idx
    return out

PFBA_FACTOR = 1.1   # as COBRApy flux_variability_analysis(pfba_factor=1.1)

def fva_tradeoff_pfba(ctx, mu, v_ref, vecs, frac=0.99, factor=PFBA_FACTOR):
    """FVA of secretion functionals with every organism's growth >= frac*mu_o
    and total absolute flux <= factor * sum|v_ref| (v_ref: parsimonious
    solution). The flux bound removes thermodynamically infeasible exchange
    cycles between organisms, which otherwise dominate the ranges."""
    n = len(ctx["rids"]); A = ctx["S"]; lb, ub = ctx["lb"], ctx["ub"]
    rows = [sp.hstack([A, -A])]; rlo = [np.zeros(A.shape[0])]; rhi = [np.zeros(A.shape[0])]
    for o, k in zip(ctx["orgs"], ctx["growth"]):
        r = np.zeros(2 * n); r[k] = 1.0; r[n + k] = -1.0
        rows.append(sp.csr_matrix(r)); rlo.append([frac * mu[o]]); rhi.append([INF])
    rows.append(sp.csr_matrix(np.ones(2 * n))); rlo.append([0.0]); rhi.append([factor * float(np.abs(v_ref).sum())])
    A2 = sp.vstack(rows).tocsc()
    lp = highspy.HighsLp(); lp.num_col_ = 2 * n; lp.num_row_ = A2.shape[0]
    lp.col_cost_ = np.zeros(2 * n)
    lp.col_lower_ = np.concatenate([np.maximum(lb, 0), np.maximum(-ub, 0)])
    lp.col_upper_ = np.concatenate([np.maximum(ub, 0), np.maximum(-lb, 0)])
    lp.row_lower_ = np.concatenate(rlo); lp.row_upper_ = np.concatenate(rhi)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A2.indptr; lp.a_matrix_.index_ = A2.indices; lp.a_matrix_.value_ = A2.data
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("simplex_strategy", 4)
    h.passModel(lp)
    out = {}; prev = np.array([], dtype=np.int32)
    for key, ks in vecs.items():
        idx = np.array([k for k, _ in ks], dtype=np.int32); val = np.array([v for _, v in ks], float)
        idx2 = np.concatenate([idx, idx + n]); val2 = np.concatenate([val, -val])
        if len(prev):
            h.changeColsCost(len(prev), prev, np.zeros(len(prev)))
        res = []
        for sgn in (1.0, -1.0):
            h.changeColsCost(len(idx2), idx2, sgn * val2); h.run()
            ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
            res.append(sgn * h.getInfo().objective_function_value if ok else np.nan)
        out[key] = (res[0], res[1]); prev = idx2
    return out
