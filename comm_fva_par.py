from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import sys, pickle, time, json, multiprocessing as mp

from comm_pipeline import *

def worker(args):
    wid, items = args
    C = Community(); C.set_medium(minimal_medium())
    res = pickle.load(open(str(RESULTS/"comm_base_res.pkl"), "rb"))
    out = {}; t0 = time.time()
    for n, (U, keys) in enumerate(items):
        out.update(run_fva(C, res, {U: keys}))
        if n % 10 == 0:
            print(f"  w{wid} {n}/{len(items)} {time.time()-t0:.0f}s", flush=True)
    return out

def main(n_workers=2):
    """Parallel FVA over the 165 contexts (worker is importable, so this also
    works with the spawn start method of macOS/Windows and inside Jupyter)."""
    t0 = time.time()
    C = Community(); C.set_medium(minimal_medium())
    res = pickle.load(open(str(RESULTS/"comm_base_res.pkl"), "rb"))
    coords = coordinates(res, C.orgs); need = needed_fva(res, coords)
    items = sorted(need.items(), key=lambda kv: -len(kv[1]))
    chunks = [items[i::n_workers] for i in range(n_workers)]
    with mp.Pool(n_workers) as pool:
        parts = pool.map(worker, list(enumerate(chunks)))
    W = {}
    for p in parts: W.update(p)
    pickle.dump(W, open(str(RESULTS/"comm_base_W.pkl"), "wb"))
    out = coherence(res, coords, W)
    print(f"T={out['T']:.2f} d={out['d']} p={out['p']:.3g} radius={out['radius']:.4f} at {out['radius_key']}", flush=True)
    for row in summarize(out, 15): print(row, flush=True)
    json.dump({"T": out["T"], "d": out["d"], "p": out["p"], "radius": out["radius"],
               "radius_key": "|".join(out["radius_key"]) if out["radius_key"] else None,
               "n_coords": len(coords), "M": sorted({k[2] for k in coords}),
               "contrib": {"|".join(k): v for k, v in out["contrib"].items()},
               "growth": {"+".join(U): {o: float(g) for o, g in r["growth"].items()} for U, r in res.items()},
               "fva_runtime_s": time.time() - t0},
              open(str(RESULTS/"comm_base.json"), "w"), indent=1)
    print("done", f"{time.time()-t0:.0f}s", flush=True)

    return out

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
