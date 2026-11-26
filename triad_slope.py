from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, json; 
from triad_lib import solve, org_exchange, pair_transfers
def point(orgs, **kw):
    com, fl, g = solve(orgs, **kw)
    bases=["glyc_e","dha_e"]; f=org_exchange(fl,orgs,bases)
    med={b:-fl.get(f"EX_{b}",0.0) for b in bases}
    x=pair_transfers(f,med,orgs,bases)
    return x[("S","H","dha_e")], -f["S"]["glyc_e"], f["S"]["dha_e"]
res={}
for name,orgs,g0,key in [("SH",("S","H"),0.0,"glyc_medium"),("DSH",("D","S","H"),5.0,"glyc_D")]:
    rows=[]
    for h in (0.25,0.5):
        lo,hi=max(g0-h,0),g0+h
        a=point(orgs,**{key:lo}); b=point(orgs,**{key:hi})
        rows.append(dict(h=h, dx_dsupply=(b[0]-a[0])/(hi-lo), dupt_dsupply=(b[1]-a[1])/(hi-lo), dsec_dsupply=(b[2]-a[2])/(hi-lo)))
    res[name]=rows; print(name,rows,flush=True)
json.dump(res,open(str(RESULTS/"triad_slope.json"),"w"),indent=1)
