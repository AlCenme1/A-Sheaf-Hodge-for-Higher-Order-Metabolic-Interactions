from paths import FIGS
from plotstyle import NAVY, RUST, TEAL, GOLD, GRAY, LIGHT, INK, savefig
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
import numpy as np, itertools

fig, ax = plt.subplots(figsize=(4.8, 4.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 9.6); ax.axis('off')

def species_dots(cx, cy, letters, color):
    n = len(letters); r = 0.32; pts = []
    for i, L in enumerate(letters):
        ang = 90 + i*(360/n)
        x = cx + r*np.cos(np.radians(ang)); y = cy + r*np.sin(np.radians(ang))
        pts.append((x, y))
        ax.add_patch(Circle((x, y), 0.14, facecolor=color, edgecolor='white', linewidth=1.0, zorder=4))
        ax.text(x, y, L, ha='center', va='center', fontsize=8, color='white', zorder=5)
    for (x1,y1),(x2,y2) in itertools.combinations(pts, 2):
        ax.plot([x1,x2], [y1,y2], color=color, lw=1.0, alpha=0.55, zorder=3)

def context_box(cx, cy, w, h, label, sublabel, color):
    box = FancyBboxPatch((cx-w/2, cy-h/2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                          facecolor=color, alpha=0.10, edgecolor=color, linewidth=1.5, zorder=2)
    ax.add_patch(box)
    ax.text(cx, cy+0.20, label, ha='center', va='center', fontsize=11, color=INK)
    ax.text(cx, cy-0.42, sublabel, ha='center', va='center', fontsize=8.7, color=GRAY)
    return cx-w/2, cx+w/2, cy-h/2, cy+h/2  # left, right, bottom, top

# --- species-dot clusters, well above the boxes ---
species_dots(2.3, 8.55, ['A','B'], NAVY)
species_dots(7.5, 8.55, ['A','B','C'], RUST)
ax.plot([2.3, 2.3], [8.05, 7.35], color=NAVY, lw=0.8, alpha=0.4, zorder=1)
ax.plot([7.5, 7.5], [8.05, 7.35], color=RUST, lw=0.8, alpha=0.4, zorder=1)

# --- context boxes ---
Ul, Ur, Ub, Ut = context_box(2.3, 6.45, 2.7, 1.7, r'$U=\{A,B\}$', r'$\mathcal{G}(U)=\mathbb{R}^{1}$', NAVY)
Vl, Vr, Vb, Vt = context_box(7.5, 6.45, 3.1, 1.7, r'$V=\{A,B,C\}$', r'$\mathcal{G}(V)=\mathbb{R}^{3}$', RUST)

# --- edge box, well BELOW both context boxes, arrows point straight down onto its top edge ---
ecx, ecy, ew, eh = 4.9, 4.55, 1.7, 0.72
eb = FancyBboxPatch((ecx-ew/2, ecy-eh/2), ew, eh, boxstyle="round,pad=0.02,rounding_size=0.06",
                     facecolor='white', edgecolor='#555555', linewidth=1.1, zorder=3)
ax.add_patch(eb)
ax.text(ecx, ecy, r'$e=\{U,V\}$', ha='center', va='center', fontsize=9.5, zorder=6)
etop = ecy + eh/2

# arrows: straight down from each box's bottom, landing on the TOP EDGE of
# the edge box (well clear of its text), never crossing through the label
ax.annotate('', xy=(ecx-0.55, etop), xytext=(2.3, Ub),
            arrowprops=dict(arrowstyle='->', color=NAVY, lw=1.4,
                             connectionstyle='arc3,rad=-0.08'))
ax.text(2.85, 5.35, r'$\rho_{U,e}$', color=NAVY, fontsize=9, ha='center')
ax.annotate('', xy=(ecx+0.55, etop), xytext=(7.5, Vb),
            arrowprops=dict(arrowstyle='->', color=RUST, lw=1.4,
                             connectionstyle='arc3,rad=0.08'))
ax.text(6.95, 5.35, r'$\rho_{V,e}$', color=RUST, fontsize=9, ha='center')

ax.text(ecx, ecy-1.05, r'$W = U\cap V = \{A,B\}$', ha='center', fontsize=9, color='#444444')
ax.text(ecx, ecy-1.48, r'coherence: $\rho_{U,e}(a_U)=\rho_{V,e}(a_V)$', ha='center', fontsize=8.6, color=INK)

ax.text(5.0, 1.55,
        r'$\mathcal{D}_0=\mathcal{U}$ (contexts)$\quad$'
        r'$\mathcal{D}_1=\{\{U,V\}:P(U\cap V)\neq\emptyset\}$ (shared pairs)',
        ha='center', fontsize=8.6, color=INK,
        bbox=dict(boxstyle='round,pad=0.35', facecolor='#F5F5F8', edgecolor='#DDDDDD', linewidth=0.8))
savefig(fig, str(FIGS/'FIG_nerve_schematic_v4.png'))
print("done")
