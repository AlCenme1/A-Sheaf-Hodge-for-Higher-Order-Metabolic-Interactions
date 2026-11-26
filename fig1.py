from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 8.5, "mathtext.fontset": "dejavuserif"})
INK, MUTED, BLUE, ORANGE = "#2b2b28", "#6b6b66", "#2a78d6", "#eb6834"
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.1, 2.75), gridspec_kw=dict(width_ratios=[1, 1.35]))
axA.set_xlim(0, 10); axA.set_ylim(0, 7.4); axB.set_xlim(0, 10); axB.set_ylim(-0.3, 6.2)
for ax in (axA, axB):
    ax.set_aspect("equal"); ax.axis("off")

def node(ax, x, y, lab, r=0.62, fc="white", ec=INK):
    ax.add_patch(Circle((x, y), r, fc=fc, ec=ec, lw=1.0, zorder=3))
    ax.text(x, y, lab, ha="center", va="center", fontsize=8.5, zorder=4, style="italic")

def arrow(ax, p, q, txt, off=(0, 0), color=INK, lw=1.2, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=10, lw=lw, color=color,
                                 shrinkA=17, shrinkB=17, connectionstyle=f"arc3,rad={rad}", zorder=2))
    mx, my = (p[0] + q[0]) / 2 + off[0], (p[1] + q[1]) / 2 + off[1]
    ax.text(mx, my, txt, ha="center", va="center", fontsize=8, color=color)

# (A) cycle on the species graph
S, H, L = (2.0, 1.4), (8.0, 1.4), (5.0, 5.8)
for p, lab, ini in ((S, "S. ruber", "$S$"), (H, "H. walsbyi", "$H$"), (L, "Halorubrum", "$L$")):
    node(axA, *p, ini)
    axA.text(p[0], p[1] - 0.95 if p[1] < 3 else p[1] + 0.95, lab, ha="center", va="center", fontsize=8, style="italic")
arrow(axA, S, L, "$2.0$", off=(-0.55, 0.15))
arrow(axA, L, H, "$1.0$", off=(0.55, 0.15))
arrow(axA, S, H, "$4.5$", off=(0, -0.35))
axA.text(5.0, 3.05, "cycle sum\n$2.0+1.0-4.5$\n$=-1.5$", ha="center", va="center", fontsize=7.5, color=MUTED, linespacing=1.3)
axA.set_title("Species graph: transitivity", fontsize=8.5, color=INK, pad=2)

# (B) context nerve: the same pair in two contexts
def ctx_box(ax, x, y, w, h, title, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.25",
                                fc="white", ec=color, lw=1.2, zorder=1))
    ax.text(x + w / 2, y + h + 0.25, title, ha="center", va="bottom", fontsize=8, color=color)
ctx_box(axB, 0.3, 0.75, 3.6, 4.45, "context $U=\\{S,H\\}$", BLUE)
ctx_box(axB, 5.6, 0.75, 4.1, 4.45, "context $V=\\{D,S,H\\}$", ORANGE)
s1, h1 = (1.1, 2.0), (3.1, 2.0)
node(axB, *s1, "$S$"); node(axB, *h1, "$H$")
arrow(axB, s1, h1, "DHA $0.8$", off=(0, 0.95), color=BLUE)
d2, s2, h2 = (7.65, 4.1), (6.5, 2.0), (8.8, 2.0)
node(axB, *d2, "$D$"); node(axB, *s2, "$S$"); node(axB, *h2, "$H$")
arrow(axB, d2, s2, "glycerol", off=(-0.75, 0.15), color=MUTED, lw=1.0)
arrow(axB, s2, h2, "DHA $1.4$", off=(0, -0.85), color=ORANGE)
axB.annotate("", xy=(5.5, 3.1), xytext=(4.0, 3.1),
             arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
axB.text(4.75, 3.45, "$W=U\\cap V$", ha="center", fontsize=8)
axB.text(4.75, 2.65, "$(\\delta_0 x)_e\\neq0$", ha="center", fontsize=8)
axB.text(5.0, 0.05, "same pair $(S,H)$, two transfers: a higher-order interaction", ha="center", fontsize=8, color=MUTED)
axB.set_title("Context nerve: context dependence", fontsize=8.5, color=INK, pad=2)
fig.tight_layout(w_pad=1.5)
fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
for ax, l in ((axA, "(A)"), (axB, "(B)")):
    tb = inv.transform(ax.get_tightbbox(r))
    fig.text(tb[0][0], tb[1][1] + 0.01, l, fontsize=10, fontweight="bold", ha="left", va="bottom")
fig.savefig(str(FIGS/"Fig1.pdf"), bbox_inches="tight")
fig.savefig(str(FIGS/"Fig1.png"), dpi=200, bbox_inches="tight")
