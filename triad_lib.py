"""Glycerol-DHA motif: joint community FBA for contexts of
D (D. salina as glycerol source), S (S. ruber iSRU780), H (H. walsbyi iHQW-TK).
"""
from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import cobra, warnings, logging, numpy as np
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import pfba, flux_variability_analysis
warnings.filterwarnings("ignore"); logging.getLogger("cobra").setLevel(logging.ERROR)

FILES = {"S": str(MODEL_S),
         "H": str(MODEL_H)}
INORG = ["h2o","h","na1","k","cl","mg2","ca2","so4","pi","nh4","fe2","fe3","mn2",
         "zn2","cu2","cobalt2","mobd","ni2","o2","co2"]
SUPP = ["gly_asn__L","lysglugly","orn","sheme"]   # essential for iSRU780

_cache = {}
def load(org):
    if org not in _cache:
        _cache[org] = cobra.io.read_sbml_model(FILES[org])
    return _cache[org]

def build_context(orgs, glyc_D=5.0, supp_rate=0.05, glyc_medium=0.0):
    """orgs subset of {'D','S','H'}. Returns cobra Model with shared pool."""
    com = Model("ctx_" + "".join(orgs))
    pool = {}
    def pool_met(base):
        if base not in pool:
            pool[base] = Metabolite(f"{base}__pool", compartment="pool")
        return pool[base]
    rxns = []
    for o in orgs:
        if o == "D":
            continue
        m = load(o).copy()
        # rename everything with suffix
        for met in m.metabolites:
            met.id = f"{met.id}__{o}"; met.compartment = f"{met.compartment}_{o}"
        for r in m.reactions:
            r.id = f"{r.id}__{o}"
        m.repair()
        for r in list(m.reactions):
            if r.id.startswith("EX_") and len(r.metabolites) == 1:
                met = list(r.metabolites)[0]
                base = met.id.rsplit("__", 1)[0]       # e.g. glyc_e
                # link: met_e(org) -> pool  (positive = secretion by organism)
                lk = Reaction(f"LINK_{base}__{o}", lower_bound=r.lower_bound, upper_bound=r.upper_bound)
                lk.add_metabolites({met: -1, pool_met(base): 1})
                rxns.append(lk)
            elif r.id.startswith("EX_"):
                continue
            else:
                rxns.append(r)
    com.add_reactions(rxns)
    # medium exchanges on the pool
    for base, pm in pool.items():
        ex = Reaction(f"EX_{base}", lower_bound=0, upper_bound=1000)
        ex.add_metabolites({pm: -1})
        com.add_reactions([ex])
    for x in INORG:
        rid = f"EX_{x}_e"
        if rid in com.reactions:
            com.reactions.get_by_id(rid).lower_bound = -1000
    for x in SUPP:
        rid = f"EX_{x}_e"
        if rid in com.reactions:
            com.reactions.get_by_id(rid).lower_bound = -supp_rate
    if glyc_medium > 0 and "EX_glyc_e" in com.reactions:
        com.reactions.get_by_id("EX_glyc_e").lower_bound = -glyc_medium
    if "D" in orgs:
        d = Reaction("SRC_glyc__D", lower_bound=glyc_D, upper_bound=glyc_D)
        d.add_metabolites({pool_met("glyc_e"): 1})
        com.add_reactions([d])
    growth = [com.reactions.get_by_id(f"Growth__{o}") for o in orgs if o != "D"]
    com.objective = {g: 1.0 for g in growth}
    return com

def org_exchange(sol_fluxes, orgs, met_bases):
    """f[o][m] = net secretion into pool by organism o (D: glycerol source)."""
    f = {}
    for o in orgs:
        f[o] = {}
        for b in met_bases:
            if o == "D":
                f[o][b] = sol_fluxes.get("SRC_glyc__D", 0.0) if b == "glyc_e" else 0.0
            else:
                f[o][b] = sol_fluxes.get(f"LINK_{b}__{o}", 0.0)
    return f

def pair_transfers(f, medium_import, orgs, met_bases, order=("D","S","H")):
    """Proportional allocation in the shared pool:
       t_{i->j,m} = s_{i,m} u_{j,m} / (sum_o s_{o,m} + medium import_m).
       Pair coordinate x_{ij,m} = t_{i->j} - t_{j->i}, i before j in `order`."""
    out = {}
    os_ = [o for o in order if o in orgs]
    for b in met_bases:
        s = {o: max(f[o][b], 0.0) for o in os_}
        u = {o: max(-f[o][b], 0.0) for o in os_}
        tot = sum(s.values()) + max(medium_import.get(b, 0.0), 0.0)
        for a in range(len(os_)):
            for c in range(a + 1, len(os_)):
                i, j = os_[a], os_[c]
                if tot <= 1e-12:
                    out[(i, j, b)] = 0.0
                else:
                    out[(i, j, b)] = (s[i] * u[j] - s[j] * u[i]) / tot
    return out

def solve(orgs, **kw):
    com = build_context(orgs, **kw)
    sol = pfba(com)
    fl = sol.fluxes.to_dict()
    growth = {o: fl.get(f"Growth__{o}", 0.0) for o in orgs if o != "D"}
    return com, fl, growth
