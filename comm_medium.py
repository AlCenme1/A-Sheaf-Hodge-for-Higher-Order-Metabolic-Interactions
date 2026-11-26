"""Organism models of the community and their growth alone on the minimal
hypersaline medium (Supplementary Table "Organism models"): reactions,
internal metabolites, exchangeable pool metabolites and growth rate (h^-1).
Glycerol 10 mmol gDW^-1 h^-1, CO2 uptake blocked, trace metals capped at 0.1.
Also checks that no organism grows without glycerol. Writes results/comm_orgs.json."""
from paths import RESULTS
import json
from comm_lib import Community, minimal_medium, fba

C = Community()
ext = set(C.ext)
rows = []
for o in C.orgs:
    R = C.org_rxns[o]
    mets = {s for r in R for s in (set(C.reactions[r]["reactants"]) | set(C.reactions[r]["products"]))}
    C.set_medium(minimal_medium()); r = fba(C.build([o]))
    mu = round(float(r["growth"][o]), 3) if r else 0.0
    C.set_medium(minimal_medium(glyc=0.0)); r0 = fba(C.build([o]))
    mu0 = float(r0["growth"][o]) if r0 else 0.0
    rows.append(dict(org=o, rxns=len(R), mets=len(mets - ext), ext=len(mets & ext), mu=mu))
    print(f"{o:18s} rxns={len(R):5d} mets={len(mets-ext):5d} ext={len(mets&ext):4d} "
          f"mu={mu:.3f}  mu(no glycerol)={mu0:.1e}", flush=True)
json.dump(rows, open(RESULTS / "comm_orgs.json", "w"), indent=1)
