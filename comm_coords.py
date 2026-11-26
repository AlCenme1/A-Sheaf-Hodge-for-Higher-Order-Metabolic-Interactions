import pickle, itertools, numpy as np
def transfers(rec, orgs, order):
    """Proportional allocation in the shared pool -> {(i,j,m): x} with i<j in `order`."""
    sec = rec["sec"]; med = rec["medium"]
    mets = {m for (_, m) in sec}
    out = {}
    os_ = sorted(orgs, key=order.index)
    for m in mets:
        s = {o: max(sec.get((o, m), 0.0), 0.0) for o in os_}
        u = {o: max(-sec.get((o, m), 0.0), 0.0) for o in os_}
        tot = sum(s.values()) + max(med.get(m, 0.0), 0.0)
        if tot <= 1e-12: continue
        for a, b in itertools.combinations(os_, 2):
            x = (s[a] * u[b] - s[b] * u[a]) / tot
            if abs(x) > 1e-9: out[(a, b, m)] = x
    return out
