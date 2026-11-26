"""Central paths. Edit here, or set the environment variables
TCBB_DATA / TCBB_RESULTS / TCBB_FIGS, instead of editing the scripts."""
import os
from pathlib import Path
ROOT = Path(__file__).resolve().parent
DATA = Path(os.environ.get("TCBB_DATA", ROOT / "data"))
RESULTS = Path(os.environ.get("TCBB_RESULTS", ROOT / "results"))
FIGS = Path(os.environ.get("TCBB_FIGS", ROOT / "figures"))
for _d in (RESULTS, FIGS):
    _d.mkdir(parents=True, exist_ok=True)

MODEL_S = DATA / "iSRU780.xml"     # Salinibacter ruber (curated CarveMe 1.6.2)
MODEL_H = DATA / "iHQW-TK.xml"     # Haloquadratum walsbyi (ISME Communications, 2026)
COMMUNITY = DATA / "community.xml" # ten-genus community (Zenodo 10.5281/zenodo.20834266)
COMM_PKL = RESULTS / "comm_parsed.pkl"   # cache of the parsed community SBML
