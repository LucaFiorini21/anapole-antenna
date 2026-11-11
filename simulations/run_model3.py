"""Model 3 - NEC-2, two generators (feeding option 2 of the report):
dipole fed at its centre + CLOSED toroidal winding with return ring, with
series capacitors that keep the current uniform, fed by its own generator.
This is the version that actually works as an anapole in a full-wave model."""
import json
import numpy as np
from scipy.optimize import minimize
from anapole.core import *
from anapole.models import *
from anapole.analysis import *
from anapole.style import *
from matplotlib.colors import LogNorm

import os
import sys
k = K0
out = {}
# "--plots-only": reuse the trimmed design and the slow results (convergence,
# sweeps, Monte Carlo, losses) stored in out/model3.json; redo only the figures.
CACHE = json.load(open("out/model3.json")) if ("--plots-only" in sys.argv and os.path.exists("out/model3.json")) else None
Pd, p0d, Zd, rd = bare_dipole()
rng = np.random.default_rng(1)


def nonuni(x, **kw):
    if min(x) <= 5:
        return 9.0
    st = model_twogen(X_turn=x[0], X_ring=x[1], **kw)
    r, V, Y = solve_ports(st, F0, {"dipole": 0.0, "coil": 1.0})
    part = np.array(r.part)
    I = r.I[(part == "coil") | (part == "return")]
    return float(np.std(I) / abs(np.mean(I)))


def tune_caps(x0=(175.0, 505.0), **kw):
    res = minimize(lambda x: nonuni(x, **kw), x0, method="Nelder-Mead", options=dict(xatol=0.05, fatol=1e-7, maxiter=400))
    return res.x, res.fun


def analyse(st, f=F0, lossy=False):
    """port basis, optimum coil/dipole current ratio alpha, tuned Result."""
    basis = port_basis(st, f, lossy)
    kk = k_of(f)
    Fd, Fc = ff(basis["dipole"][0].elements(), kk), ff(basis["coil"][0].elements(), kk)
    al = best_combination([Fd, Fc], 0)[1]
    r, V, Y = solve_ports(st, f, {"dipole": 1.0, "coil": al}, lossy)
    return dict(alpha=al, r=r, V=V, Y=Y, Fd=Fd, Fc=Fc, basis=basis)


def tuned_db(x, **kw):
    """radiated power (dB vs bare dipole) with the best current ratio, for
    capacitor reactances x[0], x[1] (ohm) and return-ring radius x[2] (mm)."""
    if x[0] < 50 or x[1] < 50:
        return 50.0
    s_ = model_twogen(X_turn=x[0], X_ring=x[1], R_ret=x[2] * 1e-3, **kw)
    b = port_basis(s_, F0)
    Fd, Fc = ff(b["dipole"][0].elements(), k), ff(b["coil"][0].elements(), k)
    al_ = best_combination([Fd, Fc], 0)[1]
    return float(db(prad_ff(Fd + al_ * Fc, k) / bare_dipole(nd=kw.get("nd", 5))[0]))


def design(x0=(175.8, 504.2, 152.71), **kw):
    """bench-style trimming: capacitors (two values) and return-ring radius
    adjusted for the deepest null; returns st, analysis, X, R_ret, non-uniformity."""
    x = np.array(x0, float)
    if CACHE and kw == dict(caps_per_turn=2, caps_ring=6):
        b_ = CACHE["baseline"]
        x = np.array([b_["X_turn"], b_["X_ring"], b_["R_ret_mm"]])
    for it in range(0 if CACHE else 2):
        sim = [x] + [x + dx for dx in ([0.5, 0, 0], [0, 5, 0], [0, 0, 0.3])]
        x = minimize(lambda y: tuned_db(y, **kw), x, method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-3, maxiter=300, initial_simplex=sim)).x
    s_ = model_twogen(X_turn=x[0], X_ring=x[1], R_ret=x[2] * 1e-3, **kw)
    return s_, analyse(s_), x[:2], x[2] * 1e-3, nonuni(x[:2], R_ret=x[2] * 1e-3, **kw)


# ---------------------------------------------------------------- 1. no capacitors
st0 = model_twogen()
r0, V0, _ = solve_ports(st0, F0, {"dipole": 0.0, "coil": 1.0})
a0 = analyse(st0)
P0 = prad_ff(ff(a0["r"].elements(), k), k)
r0e, _, _ = solve_ports(st0, F0, {"dipole": 1.0, "coil": 1.0})
out["no_caps"] = dict(equal_db=db(prad_ff(ff(r0e.elements(), k), k) / Pd), tuned_db=db(P0 / Pd), alpha=a0["alpha"],
                      Zcoil=V0["coil"])
print("no capacitors: equal currents %+.1f dB, best ratio %s -> %+.1f dB" % (out["no_caps"]["equal_db"], np.round(a0["alpha"], 3), out["no_caps"]["tuned_db"]))

# ---------------------------------------------------------------- 2. capacitor layouts
lay = []
keep = None
for cpt, cr in [(1, 6), (2, 6), (4, 6)]:
    kw = dict(caps_per_turn=cpt, caps_ring=cr)
    x, nu = tune_caps(**kw)
    st = model_twogen(X_turn=x[0], X_ring=x[1], **kw)
    re_, _, _ = solve_ports(st, F0, {"dipole": 1.0, "coil": 1.0})
    an = analyse(st)
    Pt = prad_ff(ff(an["r"].elements(), k), k)
    lay.append(dict(caps_per_turn=cpt, caps_ring=cr, n_caps=len(st.loads), X_turn=x[0], X_ring=x[1], C_turn_pF=st.C_turn * 1e12,
                    C_ring_pF=st.C_ring * 1e12, nonuni=nu, equal_db=db(prad_ff(ff(re_.elements(), k), k) / Pd),
                    alpha=an["alpha"], tuned_db=db(Pt / Pd)))
    print("caps/turn %d, ring %2d (%2d caps, %.1f / %.1f pF): non-uniformity %.2f%%, equal currents %+.1f dB, alpha %s -> %.1f dB"
          % (cpt, cr, len(st.loads), st.C_turn * 1e12, st.C_ring * 1e12, 100 * nu, lay[-1]["equal_db"], np.round(an["alpha"], 4), lay[-1]["tuned_db"]))
    if (cpt, cr) == (2, 6):
        keep = (st, an, x)
out["layouts"] = lay

out["untrimmed"] = dict(tuned_db=lay[1]["tuned_db"], alpha=lay[1]["alpha"],
                        parts=part_table(keep[1]["r"].elements(), k))
st, an, xc, RR, nu_ = design(caps_per_turn=2, caps_ring=6)
print("return ring trimmed to R_ret = %.2f mm (uniform-current value %.2f mm)" % (RR * 1e3, np.sqrt(R_T ** 2 + A_L ** 2 / 2) * 1e3))
rt = an["r"]; elt = rt.elements(); al = an["alpha"]
Ft = ff(elt, k)
out["baseline"] = dict(alpha=al, V=an["V"], tuned_db=db(prad_ff(Ft, k) / Pd), peak_db=2 * db(peak_ff(Ft) / abs(p0d)),
                       parts=part_table(elt, k), nseg=len(rt.I), X_turn=xc[0], X_ring=xc[1], R_ret_mm=RR * 1e3, nonuni=nu_,
                       C_turn_pF=st.C_turn * 1e12, C_ring_pF=st.C_ring * 1e12,
                       Zin_dipole=an["V"]["dipole"] / 1.0, Zin_coil=an["V"]["coil"] / al)
print("baseline parts:")
for n, v in out["baseline"]["parts"].items():
    print("   %-7s p0z=%8.3f mm  mz=%9.2f cm2  k2Tz=%8.3f mm  |I| %.3f..%.3f" % (n, v["p0z_mm"].real, v["mz_cm2"].real, v["k2Tz_mm"].real, v["Imin"], v["Imax"]))

# reversed winding, same port currents
stR = model_twogen(X_turn=xc[0], X_ring=xc[1], caps_per_turn=2, caps_ring=6, sense=-1, R_ret=RR)
rR, _, _ = solve_ports(stR, F0, {"dipole": 1.0, "coil": al})
out["reversed_db"] = db(prad_ff(ff(rR.elements(), k), k) / Pd)
rC, _, _ = solve_ports(st, F0, {"dipole": 0.0, "coil": al})
out["coil_alone_db"] = db(prad_ff(ff(rC.elements(), k), k) / Pd)
rD, _, _ = solve_ports(st, F0, {"dipole": 1.0, "coil": 0.0})
out["dipole_in_structure_db"] = db(prad_ff(ff(rD.elements(), k), k) / Pd)
# return ring at radius R instead of sqrt(R^2+a^2/2)
x2, _ = tune_caps(caps_per_turn=2, caps_ring=6, R_ret=R_T)
st2 = model_twogen(X_turn=x2[0], X_ring=x2[1], caps_per_turn=2, caps_ring=6, R_ret=R_T)
an2 = analyse(st2)
out["ring_at_R_db"] = db(prad_ff(ff(an2["r"].elements(), k), k) / Pd)
print("reversed %+.2f dB | coil alone %+.2f dB | return ring at R: %.1f dB" % (out["reversed_db"], out["coil_alone_db"], out["ring_at_R_db"]))

srcV = [(t, s, an["V"][g]) for g in st.ports for t, s in st.ports[g]]
if CACHE:
    for key_ in ("convergence", "sweep", "bw_khz", "cap_tolerance", "loss"):
        out[key_] = CACHE[key_]
else:
    # ---------------------------------------------------------------- 3. convergence
    conv = []
    for n_t, n_ret, nd in [(12, 36, 5), (16, 48, 5), (16, 48, 3), (24, 72, 5), (32, 96, 5)]:
        kw = dict(caps_per_turn=2, caps_ring=6, n_t=n_t, n_ret=n_ret, nd=nd)
        s_, a_, x, Rr_, _ = design((xc[0], xc[1], RR * 1e3), **kw)
        Pd_ = bare_dipole(nd=nd)[0]
        conv.append(dict(n_t=n_t, nd=nd, nseg=len(a_["r"].I), X_turn=x[0], X_ring=x[1], alpha=a_["alpha"], R_ret_mm=Rr_ * 1e3,
                         tuned_db=db(prad_ff(ff(a_["r"].elements(), k), k) / Pd_)))
        print("n_t=%2d nd=%d nseg=%3d X=(%.1f, %.1f) R_ret=%.2f mm alpha=%s tuned %.1f dB" % (n_t, nd, conv[-1]["nseg"], x[0], x[1], Rr_ * 1e3, np.round(a_["alpha"], 4), conv[-1]["tuned_db"]))
    out["convergence"] = conv

    # ---------------------------------------------------------------- 4. frequency sweeps
    fs = np.linspace(F0 - 1.5e6, F0 + 1.5e6, 301)
    srcV = [(t, s, an["V"][g]) for g in st.ports for t, s in st.ports[g]]
    fixV, fixI = [], []
    for f in fs:
        kk = k_of(f)
        Pdf = bare_dipole(f)[0]
        r = st.solve(f, sources=srcV)
        Id = r.current(*st.ports["dipole"][0])
        fixV.append(db(radiated_power(r.elements(), kk, 31, 60) / abs(Id) ** 2 / Pdf))
        r, _, _ = solve_ports(st, f, {"dipole": 1.0, "coil": al})
        fixI.append(db(radiated_power(r.elements(), kk, 31, 60) / Pdf))
    fixV, fixI = np.array(fixV), np.array(fixI)
    out["sweep"] = dict(f=fs, fixed_voltages_db=fixV, fixed_currents_db=fixI)


    def bw(y, lev):
        m = np.where(y < lev)[0]
        return float(fs[m[-1]] - fs[m[0]]) / 1e3 if len(m) else 0.0


    out["bw_khz"] = {str(l): dict(fixed_voltages=bw(fixV, l), fixed_currents=bw(fixI, l)) for l in (-20, -30, -40)}
    print("bandwidth (kHz):", out["bw_khz"])
    fig, ax = plt.subplots(figsize=(7.6, 4.1))
    ax.plot(fs / 1e6, fixI, color=COL["total"], label="generators set as current sources (I ratio fixed)")
    ax.plot(fs / 1e6, fixV, color=COL["coil"], label="generators set as voltage sources (V fixed)")
    ax.plot(fs / 1e6, 20 * np.log10(abs(1 - (fs / F0) ** 2) + 1e-6), color=COL["theory"], ls="--", label=r"report: $[1-(f/f_0)^2]^2$")
    ax.set_xlabel("frequency (MHz)"); ax.set_ylabel("radiated power vs bare dipole, same dipole current (dB)")
    ax.set_title("Bandwidth of the null (two-generator antenna, NEC-2)"); ax.legend(loc="lower right", fontsize=8.5); ax.set_ylim(-70, 25)
    save(fig, "m3_bandwidth")

    # ---------------------------------------------------------------- 5. capacitor tolerance (Monte Carlo)
    tol = {}
    base_loads = list(st.loads)
    for t_ in [0.001, 0.005, 0.01, 0.02, 0.05]:
        a_, b_ = [], []
        for _ in range(40):
            st.loads = [(tg, sg, R, L, C * (1 + rng.uniform(-t_, t_))) for tg, sg, R, L, C in base_loads]
            r = st.solve(F0, sources=srcV)
            Id = r.current(*st.ports["dipole"][0])
            a_.append(db(prad_ff(ff(r.elements(), k), k) / abs(Id) ** 2 / Pd))
            r, _, _ = solve_ports(st, F0, {"dipole": 1.0, "coil": al})
            b_.append(db(prad_ff(ff(r.elements(), k), k) / Pd))
        tol[t_] = dict(fixed_voltages=dict(median=np.median(a_), p90=np.percentile(a_, 90)),
                       fixed_currents=dict(median=np.median(b_), p90=np.percentile(b_, 90)))
        print("cap tolerance ±%.1f%%: fixed V median %.1f dB (90%% < %.1f) | fixed I median %.1f dB (90%% < %.1f)" %
              (100 * t_, np.median(a_), np.percentile(a_, 90), np.median(b_), np.percentile(b_, 90)))
    st.loads = base_loads
    out["cap_tolerance"] = tol

    # ---------------------------------------------------------------- 6. losses, efficiency, Q
    Rs = np.sqrt(np.pi * F0 * MU0 / SIGMA_CU)
    stL = model_twogen(X_turn=xc[0], X_ring=xc[1], caps_per_turn=2, caps_ring=6, R_ret=RR)
    QC = 1000.0
    stL.loads = [(tg, sg, 1 / (2 * np.pi * F0 * C) / QC, L, C) for tg, sg, R, L, C in stL.loads]
    anL = analyse(stL, lossy=True)
    rL = anL["r"]
    PradL = prad_ff(ff(rL.elements(), k), k)
    PinL = rL.p_in
    rDl, VDl, _ = solve_ports(model_loops(with_loops=False), F0, {"dipole": 1.0}, lossy=True)
    out["loss"] = dict(alpha=anL["alpha"], P_rad_W=PradL, P_in_W=PinL, efficiency=PradL / PinL, eff_db=db(PradL / PinL),
                       Z_dipole=anL["V"]["dipole"], Z_coil=anL["V"]["coil"] / anL["alpha"], tuned_db=db(PradL / Pd),
                       bare_dipole_eff=Pd / rDl.p_in, bare_dipole_Z=VDl["dipole"], Rs_mohm_sq=Rs * 1e3,
                       wire_m=stL.wire_length("coil") + stL.wire_length("return"))
    # Q of the coil circuit from the slope of its reactance (lossy), dipole port open
    df = 20e3
    z = []
    for f in (F0 - df, F0 + df):
        r, V, Y = solve_ports(stL, f, {"dipole": 0.0, "coil": 1.0}, lossy=True)
        z.append(V["coil"])
    r, V, Y = solve_ports(stL, F0, {"dipole": 0.0, "coil": 1.0}, lossy=True)
    out["loss"]["Q_coil"] = float(F0 / (2 * V["coil"].real) * (z[1].imag - z[0].imag) / (2 * df))
    out["loss"]["R_coil_ohm"] = float(V["coil"].real)
    print("lossy: P_rad/P_in = %.2e (%.1f dB), Z_dip = %s, Z_coil = %s, Q_coil = %.0f" %
          (out["loss"]["efficiency"], out["loss"]["eff_db"], np.round(out["loss"]["Z_dipole"], 3), np.round(out["loss"]["Z_coil"], 3), out["loss"]["Q_coil"]))


# ---------------------------------------------------------------- 7. figures: current, elements, patterns
part = np.array(rt.part)
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
a = ax[0]
for res, lab, c in [(a0["basis"]["coil"][0], "plain wire (no capacitors)", COL["rev"]), (an["basis"]["coil"][0], "with series capacitors", COL["coil"])]:
    pm = np.array(res.part)
    I = res.I[pm == "coil"]; s = np.arange(len(I)) / 16.0
    a.plot(s, abs(I), color=c, label=lab)
a.set_xlabel("turn number along the winding"); a.set_ylabel("|I| / I at the coil generator"); a.set_title("Current along the 12 turns")
a.legend(loc="upper right", fontsize=8.5); a.set_xticks(range(0, 13, 2))
a = ax[1]
for res, lab, c in [(a0["basis"]["coil"][0], "plain wire", COL["rev"]), (an["basis"]["coil"][0], "with capacitors", COL["coil"])]:
    pm = np.array(res.part)
    I = res.I[pm == "coil"]; s = np.arange(len(I)) / 16.0
    a.plot(s, np.degrees(np.angle(I)), color=c, label=lab)
a.set_xlabel("turn number along the winding"); a.set_ylabel("phase of I (deg)"); a.set_title("Phase along the 12 turns")
a.set_yticks([-180, -90, 0, 90, 180]); a.set_xticks(range(0, 13, 2))
save(fig, "m3_coil_current")

# contributions of each element to the far field on the equator (bridge chart)
th90, ph15 = np.array([np.pi / 2]), np.array([np.radians(15)])
grp = [("dipole + plates", elt.sel("dipole", "plates"), COL["dipole"]), ("toroidal winding", elt.sel("coil", "leads"), COL["coil"]),
       ("return ring", elt.sel("return"), COL["ret"])]
contrib = []
for lab, e, c in grp:
    Nt, Np = far_field(e, th90, ph15, k)
    contrib.append((lab, complex(-Nt[0]) * 1e3, complex(Np[0]) * 1e3, c))
Ntot = far_field(elt, th90, ph15, k)
out["baseline"]["contrib_equator_mm"] = {l: dict(theta=t, phi=p) for l, t, p, c in contrib}
out["baseline"]["contrib_equator_mm"]["sum"] = dict(theta=complex(-Ntot[0][0]) * 1e3, phi=complex(Ntot[1][0]) * 1e3)
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
for a, key, ttl, part_ in [(ax[0], 1, r"$\theta$-polarized field (in phase with the dipole current)", "real"),
                           (ax[1], 2, r"$\varphi$-polarized field (in quadrature)", "imag")]:
    run_ = 0.0
    for i, cc in enumerate(contrib):
        v = getattr(cc[key], part_)
        a.bar(i, v, bottom=run_, color=cc[3], width=0.55)
        a.text(i, max(run_, run_ + v) + (1.0 if key == 1 else 4), "%+.2f" % v, ha="center", va="bottom", fontsize=9, color=INK2)
        a.plot([i + 0.275, i + 0.725], [run_ + v] * 2, color=MUTED, lw=0.8)
        run_ += v
    a.bar(len(contrib), run_, color=COL["total"], width=0.55)
    a.text(len(contrib), run_ + (1.0 if key == 1 else 4), "%+.3f" % run_, ha="center", va="bottom", fontsize=9, color=INK)
    a.set_xticks(range(len(contrib) + 1)); a.set_xticklabels(["dipole\n+ plates", "toroidal\nwinding", "return\nring", "sum"], fontsize=9)
    a.axhline(0, color=AXIS, lw=1); a.grid(axis="x", visible=False); a.set_title(ttl)
ax[0].set_ylabel("far field on the equator, as effective length (mm)")
ax[0].set_ylim(-12, 36); ax[1].set_ylim(-40, 160)
save(fig, "m3_element_contributions")

# far-field cuts of each element (NEC currents)
th = np.linspace(0, np.pi, 181)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for lab, e, c in grp + [("sum of all parts", elt, COL["total"])]:
    Nt, Np = far_field(e, th, 0 * th + np.radians(15), k)
    ax[0].plot(np.degrees(th), -Nt.real / abs(p0d), color=c, label=lab, lw=2.4 if c == COL["total"] else 2)
    ax[1].plot(np.degrees(th), Np.imag / abs(p0d), color=c, label=lab, lw=2.4 if c == COL["total"] else 2)
ax[0].set_title(r"$\theta$-polarized far field of each element"); ax[1].set_title(r"$\varphi$-polarized far field of each element")
ax[0].set_ylabel(r"signed field / field of the bare dipole at $\theta$ = 90°")
for a in ax:
    a.set_xlabel(r"$\theta$ (deg, from the dipole axis)"); a.set_xticks(range(0, 181, 30)); a.axhline(0, color=AXIS, lw=1)
ax[1].legend(loc="upper right", fontsize=8.5); ax[1].set_ylim(-5, 7.5)
save(fig, "m3_element_fields")

# 3D patterns
T, P = np.meshgrid(np.linspace(0, np.pi, 61), np.linspace(0, 2 * np.pi, 121), indexing="ij")
pat = [("dipole alone", rd.elements(), COL["dipole"]), ("coil alone (same currents)", rC.elements(), COL["coil"]),
       ("dipole + coil (anapole)", elt, COL["total"]), ("dipole + reversed coil", rR.elements(), COL["rev"])]
fig = plt.figure(figsize=(13, 3.9))
for i, (name, e, c) in enumerate(pat):
    Nt, Np = far_field(e, T, P, k)
    A = np.sqrt(abs(Nt) ** 2 + abs(Np) ** 2) / abs(p0d)
    a = fig.add_subplot(1, 4, i + 1, projection="3d")
    a.plot_surface(A * np.sin(T) * np.cos(P), A * np.sin(T) * np.sin(P), A * np.cos(T), color=c, alpha=0.85, linewidth=0.15,
                   edgecolor=SURF, rstride=2, cstride=4, shade=True)
    a.plot([0, 0], [0, 0], [-1.6, 1.6], color=MUTED, lw=0.8)
    a.set_xlim(-2, 2); a.set_ylim(-2, 2); a.set_zlim(-1.6, 1.6); a.set_box_aspect((1, 1, 0.8)); a.set_axis_off(); a.view_init(22, -55)
    a.set_title("%s\npeak field %.4f × dipole (%+.1f dB)" % (name, A.max(), 20 * np.log10(A.max())), fontsize=8.5, loc="center", fontweight="normal")
fig.subplots_adjust(wspace=-0.05)
save(fig, "m3_patterns_3d")

# residual map
Tg, Pg = np.meshgrid(np.linspace(0, np.pi, 91), np.linspace(0, 2 * np.pi, 181), indexing="ij")
Nt, Np = far_field(elt, Tg, Pg, k)
D = 20 * np.log10(np.sqrt(abs(Nt) ** 2 + abs(Np) ** 2) / abs(p0d) + 1e-12)
fig, a = plt.subplots(figsize=(7.2, 3.8))
im = a.pcolormesh(np.degrees(Pg), np.degrees(Tg), D, cmap="Blues", vmin=-80, vmax=-30, shading="auto")
a.set_xlabel(r"$\varphi$ (deg)"); a.set_ylabel(r"$\theta$ (deg)"); a.invert_yaxis(); a.grid(False); a.set_xticks(range(0, 361, 60))
a.set_title("Residual far field of the tuned antenna — max %.1f dB" % D.max())
fig.colorbar(im, ax=a, pad=0.02).set_label("dB (peak of the bare dipole = 0)")
save(fig, "m3_residual_map")
out["baseline"]["residual_max_db"] = float(D.max())

# ---------------------------------------------------------------- 8. near field from the NEC currents
ext, n = 0.6, 181
x = np.linspace(-ext, ext, n); X, Zz = np.meshgrid(x, x, indexing="xy")
pts = {"xz": np.stack([X, 0 * X + 0.004, Zz], -1), "xy": np.stack([X, Zz, 0 * X + 0.012], -1)}
sets = [("dipole alone", elt.sel("dipole", "plates")), ("winding + return ring alone", elt.sel("coil", "leads", "return")), ("all parts (anapole)", elt)]
tt = np.linspace(0, 2 * np.pi, 100)
for comp, idx, unit, vmin, vmax, pl, ttl in [("E", 0, "V/m", 1e-1, 1e5, "xz", "Electric field |E|, meridian plane (x–z), NEC-2 currents"),
                                             ("H", 1, "A/m", 1e-4, 1e2, "xy", "Magnetic field |H|, equatorial plane (x–y), NEC-2 currents")]:
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
    for a, (lab, e) in zip(ax, sets):
        EH = near_field(e, pts[pl].reshape(-1, 3), k)
        A = np.linalg.norm(EH[idx], axis=1).reshape(n, n)
        im = a.pcolormesh(x, x, np.clip(A, vmin, vmax), norm=LogNorm(vmin, vmax), cmap="magma", shading="auto")
        a.plot((R_T + A_L) * np.cos(tt), (R_T + A_L) * np.sin(tt), color="white", lw=0.6, ls=":")
        a.set_aspect("equal"); a.grid(False); a.set_title(lab); a.set_xlabel("x (m)")
    ax[0].set_ylabel("z (m)" if pl == "xz" else "y (m)")
    fig.colorbar(im, ax=ax, pad=0.015, shrink=0.9).set_label("|%s| (%s) for 1 A in the dipole" % (comp, unit))
    fig.suptitle(ttl + " — dotted: sphere r = R + a", x=0.07, ha="left", fontsize=11, color=INK)
    save(fig, "m3_near_%s_%s" % (comp, pl))

rr = np.geomspace(0.03, 30, 120)
u = np.array([np.cos(np.radians(15)), np.sin(np.radians(15)), 0.0])
fig, a = plt.subplots(figsize=(7.6, 4.2))
dec = {}
for (lab, e), c in zip(sets, [COL["dipole"], COL["coil"], COL["total"]]):
    E, H = near_field(e, rr[:, None] * u[None], k)
    dec[lab] = np.linalg.norm(E, axis=1)
    a.loglog(rr, dec[lab], color=c, label=lab)
a.axvline(R_T + A_L, color=AXIS, lw=1); a.text((R_T + A_L) * 1.1, 2e4, "r = R + a", color=INK2, fontsize=9)
a.axvline(LAM0 / (2 * np.pi), color=AXIS, lw=1, ls=":"); a.text(LAM0 / (2 * np.pi) * 1.08, 2e3, r"r = $\lambda/2\pi$", color=INK2, fontsize=9)
a.set_xlabel("distance from the centre r (m), on the equator"); a.set_ylabel("|E| (V/m) for 1 A in the dipole")
a.set_title("NEC-2: the field stays near the antenna"); a.legend(loc="lower left"); a.set_ylim(1e-6, 1e6)
save(fig, "m3_decay")
key = list(dec)
out["decay"] = dict(r=rr, **{kk_: v for kk_, v in dec.items()})
for r_ in (0.3, 0.5, 1.0, 3.0, 30.0):
    i = np.argmin(abs(rr - r_))
    out.setdefault("suppression_vs_r_db", {})[str(r_)] = float(20 * np.log10(dec[key[2]][i] / dec[key[0]][i]))
print("field suppression vs distance (dB):", out["suppression_vs_r_db"])

# ---------------------------------------------------------------- decks
st.sources = srcV
open("nec/model3_two_generators_tuned.nec", "w").write(
    st.deck(F0, "Anapole - dipole + closed toroidal winding with series capacitors, two generators, tuned"))
st0.sources = [(t, s, V0g) for g, V0g in a0["V"].items() for t, s in st0.ports[g]]
open("nec/model3_two_generators_no_capacitors.nec", "w").write(st0.deck(F0, "Two generators, plain winding without capacitors"))
json.dump(jsonable(out), open("out/model3.json", "w"), indent=1)
