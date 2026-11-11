"""Model 1 - NEC-2, dipole with end plates + 12 independent loops.
The geometry of the report with REAL current distributions: what changes
when the loop current is not perfectly uniform, and how deep the null is
once the feed currents are tuned."""
import json
import numpy as np
from anapole.core import *
from anapole.models import *
from anapole.analysis import *
from anapole.style import *

k = K0
out = {}
Pd, p0d, Zd, rd = bare_dipole()
out["bare_dipole"] = dict(Rrad_mohm=2 * Pd * 1e3, h_eff_mm=p0d.real * 1e3, Zin=Zd, P_in_over_P_rad=rd.p_in / Pd,
                          I_segments=np.abs(rd.I[np.array(rd.part) == "dipole"]))
print("bare dipole with plates: Rrad = %.2f mOhm, h_eff = %.2f mm, Zin = %.3f %+.0fj ohm" %
      (2 * Pd * 1e3, p0d.real * 1e3, Zd.real, Zd.imag))

# ---- loop alone: current distribution for the feed variants -------------
loopI = {}
for feed in ["outer", "inner", "top", "dual", "quad"]:
    st = model_loops(with_dipole=False, feed=feed, n_l=32)
    r, V, Y = solve_ports(st, F0, {g: 1.0 for g in st.ports})
    loopI[feed] = r.I[:32]
    out.setdefault("loop_alone", {})[feed] = dict(Z=V["loops"], Imin=abs(r.I[:32]).min(), Imax=abs(r.I[:32]).max(),
                                                   parts=part_table(r.elements(), k))

# ---- anapole with equal feed currents, then tuned ------------------------
cases = {}
for feed in ["outer", "inner", "top", "dual", "quad"]:
    st = model_loops(feed=feed)
    groups = list(st.ports)
    r, V, Y = solve_ports(st, F0, {g: 1.0 for g in groups})
    el = r.elements()
    F = ff(el, k)
    basis = port_basis(st, F0)
    Fs = [ff(basis[g][0].elements(), k) for g in groups]
    w2 = best_combination(Fs, 0)                               # every port group free (exotic)
    P2 = prad_ff(sum(wi * Fi for wi, Fi in zip(w2, Fs)), k)
    al = best_combination([Fs[0], sum(Fs[1:])], 0)[1]          # ONE common factor for all loop feeds
    w = np.array([1.0] + [al] * (len(groups) - 1))
    rt, Vt, _ = solve_ports(st, F0, dict(zip(groups, w)))
    elt = rt.elements(); Ft = ff(elt, k)
    cases[feed] = dict(equal=dict(P_db=db(prad_ff(F, k) / Pd), peak_db=2 * db(peak_ff(F) / abs(p0d)), parts=part_table(el, k)),
                       tuned=dict(P_db=db(prad_ff(Ft, k) / Pd), peak_db=2 * db(peak_ff(Ft) / abs(p0d)), w=w, V=Vt,
                                  parts=part_table(elt, k)),
                       free=dict(P_db=db(P2 / Pd), w=w2))
    print("feed=%-5s equal currents: P/Pdip = %6.1f dB | tuned currents %s -> %6.1f dB (peak %.1f dB)" %
          (feed, cases[feed]["equal"]["P_db"], np.round(w[1:], 3), cases[feed]["tuned"]["P_db"], cases[feed]["tuned"]["peak_db"]))
    if feed == "dual":
        keep = (st, rt, el, elt, w, basis, groups)
out["cases"] = cases

# reversed winding (dual feed, equal currents)
st = model_loops(feed="dual", sense=-1)
r, V, Y = solve_ports(st, F0, {g: 1.0 for g in st.ports})
out["reversed_db"] = db(prad_ff(ff(r.elements(), k), k) / Pd)
print("reversed winding, equal currents: %+.2f dB" % out["reversed_db"])

# ---- convergence ---------------------------------------------------------
conv = []
for n_l in [8, 12, 16, 24, 32, 48]:
    for nd in [3, 5, 7]:
        if nd != 5 and n_l != 16:
            continue
        kw = dict(nd=nd)
        Pd_, p0_, _, _ = bare_dipole(nd=nd)
        st = model_loops(feed="dual", n_l=n_l, **kw); groups = list(st.ports)
        r, V, Y = solve_ports(st, F0, {g: 1.0 for g in groups})
        Pe = prad_ff(ff(r.elements(), k), k)
        basis = port_basis(st, F0)
        Fs = [ff(basis[g][0].elements(), k) for g in groups]
        al = best_combination([Fs[0], sum(Fs[1:])], 0)[1]
        w = np.array([1.0, al, al])
        Pt = prad_ff(sum(wi * Fi for wi, Fi in zip(w, Fs)), k)
        m = moments(basis["loops"][0].elements(), k)
        conv.append(dict(n_l=n_l, nd=nd, nseg=len(r.I), equal_db=db(Pe / Pd_), tuned_db=db(Pt / Pd_), w=w, h_eff_mm=p0_.real * 1e3))
        print("n_l=%2d nd=%d nseg=%4d equal %.1f dB tuned %.1f dB  w=%s" % (n_l, nd, len(r.I), conv[-1]["equal_db"], conv[-1]["tuned_db"], np.round(w[1:], 4)))
out["convergence"] = conv

# ---- figures ---------------------------------------------------------------
# (a) loop current distribution vs feed position
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw=dict(width_ratios=[1.2, 1], wspace=0.28))
psi = np.arange(32) * 360 / 32
a = ax[0]
for feed, c, lab in [("outer", COL["coil"], "fed on the outer side (ψ = 0°)"), ("inner", COL["dipole"], "fed on the inner side (ψ = 180°)"),
                     ("top", COL["ret"], "fed at the top (ψ = 90°)"), ("dual", COL["total"], "two feeds: outer + inner")]:
    a.plot(psi, abs(loopI[feed]), color=c, label=lab, marker="o", ms=3)
a.set_xlabel("position around the loop ψ (deg); 0° = outer side, 180° = side facing the dipole")
a.set_ylabel("|I| / I at the feed"); a.set_xticks(range(0, 361, 45)); a.set_title("One loop is not a uniform current")
a.legend(loc="upper left", fontsize=8, ncol=2); a.set_ylim(0.99, 1.1)
a = ax[1]
feeds = ["outer", "inner", "top", "quad", "dual"]
ye = [cases[f_]["equal"]["P_db"] for f_ in feeds]; yt = [cases[f_]["tuned"]["P_db"] for f_ in feeds]
y = np.arange(len(feeds))
a.barh(y - 0.19, ye, 0.34, color=COL["coil"], label="equal feed currents (report)")
a.barh(y + 0.19, yt, 0.34, color=COL["total"], label="feed currents tuned")
for yi, v in zip(y - 0.19, ye):
    a.text(v + (1.2 if v >= 0 else -1.2), yi, "%+.1f" % v, va="center", ha="left" if v >= 0 else "right", fontsize=8.5, color=INK2)
for yi, v in zip(y + 0.19, yt):
    a.text(v - 1.2, yi, "%.1f" % v, va="center", ha="right", fontsize=8.5, color=INK2)
a.set_yticks(y); a.set_yticklabels(["outer", "inner", "top", "4 feeds", "outer\n+ inner"]); a.invert_yaxis()
a.axvline(0, color=AXIS, lw=1); a.set_xlim(-75, 12); a.grid(axis="y", visible=False)
a.set_xlabel("radiated power vs bare dipole (dB)"); a.set_title("Feed point of the loops vs depth of the null")
a.legend(loc="upper center", fontsize=8.5, bbox_to_anchor=(0.5, -0.17), ncol=2)
save(fig, "m1_feed_position")

# (b) element fields from the NEC currents (dual feed, tuned)
st, rt, el, elt, w, basis, groups = keep
th = np.linspace(0, np.pi, 181)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for a, (e, ttl) in zip(ax, [(el, "NEC-2, equal feed currents"), (elt, "NEC-2, feed currents tuned")]):
    for n, lab in [("dipole", "dipole + plates"), ("coil", "12 loops"), ("total", "sum")]:
        ee = e if n == "total" else (e.sel("dipole", "plates") if n == "dipole" else e.sel("coil"))
        Nt, Np = far_field(ee, th, 0 * th + np.radians(15), k)
        a.plot(np.degrees(th), -Nt.real / abs(p0d), color=COL[n], label=lab)
    a.set_title(ttl); a.set_xlabel(r"$\theta$ (deg)"); a.set_xticks(range(0, 181, 30)); a.axhline(0, color=AXIS, lw=1)
ax[0].set_ylabel(r"signed far field  $-$Re $N_\theta$ / ($I\,h_{\rm eff}$)"); ax[0].legend(loc="upper right")
save(fig, "m1_element_fields")
json.dump(jsonable(out), open("out/model1.json", "w"), indent=1)

# decks for 4nec2
st.sources = rt.src
open("nec/model1_loops_dualfeed_tuned.nec", "w").write(st.deck(F0, "Anapole - dipole + 12 independent loops, dual feed, tuned currents"))
s0 = model_loops(with_loops=False); s0.source(*s0.ports["dipole"][0], 1.0)
open("nec/reference_dipole_with_plates.nec", "w").write(s0.deck(F0, "Reference: 29.8 mm dipole with 30 mm end plates"))
