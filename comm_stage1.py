from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, pickle, itertools; 
from comm_pipeline import *
C = Community(); C.set_medium(minimal_medium())
contexts = [c for k in (2, 3) for c in itertools.combinations(C.orgs, k)]
res = solve_contexts(C, contexts, log=True)
pickle.dump(res, open(str(RESULTS/"comm_base_res.pkl"), "wb")); print("stage1 done", flush=True)
