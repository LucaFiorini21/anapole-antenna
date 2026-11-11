"""3-D drawings of the wire models."""
import numpy as np
from anapole.core import *
from anapole.models import *
from anapole.style import *

mods = [("Model 1 — 12 independent loops", model_loops(feed="dual")),
        ("Model 2 — series, one generator", model_series()),
        ("Model 3 — two generators, capacitors", model_twogen(caps_per_turn=2, X_turn=175.8, caps_ring=6, X_ring=504.2))]
cmap = dict(dipole=COL["dipole"], plates=COL["dipole"], coil=COL["coil"], leads=COL["leads"])
cmap["return"] = COL["ret"]
fig = plt.figure(figsize=(13.5, 3.9))
for i, (ttl, st) in enumerate(mods):
    a = fig.add_subplot(1, 3, i + 1, projection="3d")
    for w in st.wires:
        p, q = w["p1"], w["p2"]
        a.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], color=cmap[w["part"]], lw=1.6 if w["part"] != "plates" else 1.0)
    tags = {w["tag"]: w for w in st.wires}
    for t, sg, *_ in st.loads:
        c = 0.5 * (tags[t]["p1"] + tags[t]["p2"]); a.scatter(*c, color=INK, s=9, depthshade=False)
    for g in st.ports.values():
        for t, sg in g:
            w = tags[t]; c = w["p1"] + (sg - 0.5) / w["nseg"] * (w["p2"] - w["p1"])
            a.scatter(*c, color=COL["rev"], s=26, marker="D", depthshade=False)
    a.set_xlim(-0.2, 0.2); a.set_ylim(-0.2, 0.2); a.set_zlim(-0.1, 0.1); a.set_box_aspect((1, 1, 0.5), zoom=1.45); a.set_axis_off()
    a.view_init(28, -62); a.set_title(ttl, fontsize=10, loc="center")
h = [plt.Line2D([], [], color=COL["dipole"], lw=2, label="dipole with end plates"), plt.Line2D([], [], color=COL["coil"], lw=2, label="toroidal coil"),
     plt.Line2D([], [], color=COL["ret"], lw=2, label="return ring"), plt.Line2D([], [], color=COL["leads"], lw=2, label="leads"),
     plt.Line2D([], [], color=COL["rev"], marker="D", ls="", label="generator"), plt.Line2D([], [], color=INK, marker="o", ms=3, ls="", label="series capacitor")]
fig.legend(handles=h, loc="lower center", ncol=6, fontsize=9)
fig.subplots_adjust(wspace=-0.05, bottom=0.08)
save(fig, "geometry")
