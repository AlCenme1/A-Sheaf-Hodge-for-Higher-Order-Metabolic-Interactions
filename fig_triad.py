from paths import DATA, RESULTS, FIGS, MODEL_S, MODEL_H, COMMUNITY, COMM_PKL
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":8,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.edgecolor":"#6b6b66","axes.linewidth":0.6,"xtick.color":"#4a4a46","ytick.color":"#4a4a46"})
b=json.load(open(str(RESULTS/"triad_base.json")))["contrib"]
coords=[("D|S|glyc_e","DS","D–S\nglycerol"),("D|H|glyc_e","DH","D–H\nglycerol"),("S|H|dha_e","SH","S–H\nDHA")]
C_PAIR,C_TRI="#2a78d6","#eb6834"
fig,(axA,axB)=plt.subplots(1,2,figsize=(7.0,2.6),gridspec_kw=dict(width_ratios=[1.15,1]))
x=np.arange(len(coords))
for k,(key,pc,lab) in enumerate(coords):
    v=b[key]["values"]; s=b[key]["sigma"]
    axA.errorbar(k-0.12,v[pc],yerr=s[pc],fmt="o",ms=5,color=C_PAIR,capsize=2,lw=1,label="Pair context" if k==0 else None)
    axA.errorbar(k+0.12,v["DSH"],yerr=s["DSH"],fmt="s",ms=5,color=C_TRI,capsize=2,lw=1,label="Triad {D,S,H}" if k==0 else None)
axA.set_xticks(x); axA.set_xticklabels([c[2] for c in coords])
axA.set_ylabel("Net transfer $x_{ij,m}$\n(mmol gDW$^{-1}$ h$^{-1}$)")
axA.axhline(0,color="#c3c2b7",lw=0.6,zorder=0); axA.grid(axis="y",color="#e8e7e2",lw=0.5); axA.set_axisbelow(True)
axA.legend(frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.22),ncol=2)
top=sorted(b.items(),key=lambda kv:-kv[1]["T"])[:6]
labs=[k.replace("|","–",1).replace("|"," ").replace("_e","").replace("__L","").replace("glyc","glycerol").replace("dha","DHA") for k,_ in top]
vals=[v["T"] for _,v in top]
y=np.arange(len(top))[::-1]
axB.barh(y,vals,height=0.6,color="#2a78d6")
axB.set_xscale("log"); axB.set_yticks(y); axB.set_yticklabels(labs)
axB.set_xlabel("Contribution to T (log scale)"); axB.grid(axis="x",color="#e8e7e2",lw=0.5); axB.set_axisbelow(True)
for yy,vv in zip(y,vals): axB.text(vv*1.15,yy,f"{vv:.2f}",va="center",fontsize=7,color="#4a4a46")
axB.set_xlim(min(vals)/2,max(vals)*8)
fig.tight_layout(w_pad=2.5)
fig.canvas.draw()
r=fig.canvas.get_renderer(); inv=fig.transFigure.inverted()
for ax,l in ((axA,"(A)"),(axB,"(B)")):
    tb=inv.transform(ax.get_tightbbox(r)); bb=ax.get_position()
    fig.text(tb[0][0],bb.y1+0.05,l,fontsize=10,fontweight="bold",ha="left",va="bottom")
fig.savefig(str(FIGS/"FIG_triad.pdf"),bbox_inches="tight"); fig.savefig(str(FIGS/"FIG_triad.png"),dpi=200,bbox_inches="tight")
