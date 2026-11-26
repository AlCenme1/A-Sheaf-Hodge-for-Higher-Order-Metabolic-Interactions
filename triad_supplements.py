"""Essential supplements of iSRU780 (Section VII-A). Starting from a complete
medium with glycerol 10 and organic uptakes capped at 1, remove the organic
components one at a time, in alphabetical order, while growth persists.
Expected result: Gly-Asn, Lys-Glu-Gly, ornithine and siroheme.
Writes results/sr_supp.json."""
from paths import RESULTS, MODEL_S
import cobra, warnings, logging, json
warnings.filterwarnings("ignore"); logging.getLogger("cobra").setLevel(logging.ERROR)
from triad_lib import INORG

m = cobra.io.read_sbml_model(str(MODEL_S))
ex = [r.id for r in m.exchanges]
base = {f"EX_{x}_e": 1000 for x in INORG if f"EX_{x}_e" in ex}
med = {r: 1000 for r in ex}; med["EX_glyc_e"] = 10
for r in ex:
    if r not in base and r != "EX_glyc_e":
        med[r] = 1                      # organics cannot act as the carbon source
m.medium = med
print("complete medium, organics capped at 1:", round(m.slim_optimize(), 4))
keep = dict(med)
for r in sorted(r for r in ex if r not in base and r != "EX_glyc_e"):
    trial = dict(keep); trial.pop(r, None)
    m.medium = trial
    if m.slim_optimize(error_value=0) > 1e-3:
        keep = trial
supp = [r for r in keep if r not in base and r != "EX_glyc_e"]
m.medium = keep
print("essential supplements:", supp, "growth", round(m.slim_optimize(), 4))
json.dump(supp, open(RESULTS / "sr_supp.json", "w"))
