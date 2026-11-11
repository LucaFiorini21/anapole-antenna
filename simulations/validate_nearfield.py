"""Cross-check: near field from my Hertzian-dipole sum of the NEC segment
currents vs NEC-2's own near-field routine (NE/NH cards), series model."""
import json
import numpy as np
from PyNEC import nec_context
from anapole.core import *
from anapole.models import *
st = model_twogen(X_turn=175.8, X_ring=504.2, R_ret=0.15271, caps_per_turn=2, caps_ring=6)
r, V, Y = solve_ports(st, F0, {"dipole": 1.0, "coil": 0.7655})
ctx = nec_context(); g = ctx.get_geometry()
for w in st.wires:
    a, b = w["p1"], w["p2"]; g.wire(w["tag"], w["nseg"], *a, *b, w["rad"], 1.0, 1.0)
ctx.geometry_complete(0); ctx.set_extended_thin_wire_kernel(True)
for t, sg, R, L, C in st.loads:
    ctx.ld_card(0, t, sg, sg, R, L, C)
ctx.fr_card(0, 1, F0 / 1e6, 0)
for t, sg, Vv in r.src:
    ctx.ex_card(0, t, sg, 0, Vv.real, Vv.imag, 0, 0, 0, 0)
xs = np.array([0.05, 0.25, 0.4, 0.8, 2.0, 6.0])
rows = []
for x in xs:
    ctx.ne_card(0, 1, 1, 1, x, 0.03, 0.02, 0, 0, 0)
    nf = ctx.get_near_field_pattern(len(rows) * 2)
    En = np.array([nf.get_field_x()[0], nf.get_field_y()[0], nf.get_field_z()[0]])
    ctx.nh_card(0, 1, 1, 1, x, 0.03, 0.02, 0, 0, 0)
    nf = ctx.get_near_field_pattern(len(rows) * 2 + 1)
    Hn = np.array([nf.get_field_x()[0], nf.get_field_y()[0], nf.get_field_z()[0]])
    E, H = near_field(r.elements(), np.array([[x, 0.03, 0.02]]), r.k)
    rows.append(dict(x=x, E_nec=float(np.linalg.norm(En)), E_mine=float(np.linalg.norm(E)), H_nec=float(np.linalg.norm(Hn)), H_mine=float(np.linalg.norm(H)),
                     dE=float(np.linalg.norm(E[0] - En) / np.linalg.norm(En)), dH=float(np.linalg.norm(H[0] - Hn) / np.linalg.norm(Hn))))
    print("x=%.2f m  |E| NEC %.4e mine %.4e (vector diff %.2f%%)   |H| NEC %.4e mine %.4e (%.2f%%)" %
          (x, rows[-1]["E_nec"], rows[-1]["E_mine"], 100 * rows[-1]["dE"], rows[-1]["H_nec"], rows[-1]["H_mine"], 100 * rows[-1]["dH"]))
json.dump(rows, open("out/validate_nearfield.json", "w"), indent=1)
