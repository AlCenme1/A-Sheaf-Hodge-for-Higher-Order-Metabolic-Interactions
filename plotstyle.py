"""Shared professional plotting style for the paper's figures.
Uses Latin Modern (CM-like) fonts to match IEEEtran's typography,
a restrained academic palette, and no in-figure titles (captions live
in LaTeX \\caption{} instead)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAVY   = '#1B2A4A'
RUST   = '#B5533C'
TEAL   = '#2E7D6B'
GOLD   = '#B08D2B'
GRAY   = '#8C8C8C'
LIGHT  = '#E7E7EC'
INK    = '#1A1A1A'

PALETTE = [NAVY, RUST, TEAL, GOLD]

plt.rcParams.update({
    'font.family': 'Latin Modern Roman',
    'mathtext.fontset': 'cm',
    'font.size': 10,
    'axes.labelsize': 10.5,
    'axes.linewidth': 0.8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8.7,
    'legend.frameon': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'axes.edgecolor': '#444444',
    'xtick.color': '#333333',
    'ytick.color': '#333333',
    'text.color': INK,
    'axes.labelcolor': INK,
})

def style_ax(ax, grid_axis='y'):
    ax.set_facecolor('white')
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=LIGHT, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color('#444444')
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(axis='both', length=3, width=0.7)

def savefig(fig, path, dpi=300):
    fig.savefig(path, dpi=dpi, bbox_inches='tight', pad_inches=0.03)
