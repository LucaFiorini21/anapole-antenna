"""Model 0 - ideal (prescribed) currents. Reproduces the hypotheses of the
report with exact field integrals (no multipole expansion) and shows the field
of every element and of their sum, in the far and in the near zone."""
import json
import numpy as np
from anapole.core import *
from anapole import ideal
from anapole.style import *
from matplotlib.colors import LogNorm

k = K0
out = {}
R_RET = np.sqrt(R_T ** 2 + A_L ** 2 / 2)      # return-ring radius that cancels m_z of the helix
out["R_ret_exact_mm"] = R_RET * 1e3


def peak_db(e, kk=k, ref=H_D):
    T, P, W = sphere_grid(61, 120)
    Nt, Np = far_field(e, T, P, kk)
    return 20 * np.log10(np.sqrt(abs(Nt) ** 2 + abs(Np) ** 2).max() / ref)


# ---------------------------------------------------------------- table 1
rows = []
cases = [
    ("dipole alone", ideal.dipole()),
    ("coil alone: 12 separate loops", ideal.loops()),
    ("coil alone: helix, no return wire", ideal.helix()),
    ("return ring alone (radius R)", ideal.ring()),
    ("helix + return ring at R", Elements.concat([ideal.helix(), ideal.ring()])),
    ("helix + return ring at sqrt(R^2+a^2/2)", Elements.concat([ideal.helix(), ideal.ring(R=R_RET)])),
    ("ANAPOLE: dipole + 12 loops", ideal.anapole_loops()),
    ("ANAPOLE, h tuned to 29.46 mm", ideal.anapole_loops(h=0.02946)),
    ("dipole + 12 loops, winding reversed", ideal.anapole_loops(sense=-1)),
    ("dipole + helix, NO return wire", ideal.anapole_helix(with_return=False)),
    ("dipole + helix + return ring at R", ideal.anapole_helix()),
    ("dipole + helix + return ring at sqrt(R^2+a^2/2)",
     Elements.concat([ideal.dipole(), ideal.helix(), ideal.ring(R=R_RET)])),
]
P_dip = radiated_power(ideal.dipole(), k)
for name, e in cases:
    m = moments(e, k)
    P = radiated_power(e, k)
    rows.append(dict(case=name, p0z_mm=m["p0"][2].real * 1e3, mz_cm2=m["m"][2].real * 1e4,
                     k2Tz_mm=k * k * m["T"][2].real * 1e3, Rrad_mohm=2 * P * 1e3,
                     P_rel_db=10 * np.log10(P / P_dip), peak_db=peak_db(e)))
    print("%-48s p0z=%7.3f mm  mz=%8.3f cm2  k2Tz=%8.3f mm  Rrad=%10.5f mOhm  P/Pdip=%7.2f dB  peak=%6.1f dB"
          % (name, rows[-1]["p0z_mm"], rows[-1]["mz_cm2"], rows[-1]["k2Tz_mm"], rows[-1]["Rrad_mohm"],
             rows[-1]["P_rel_db"], rows[-1]["peak_db"]))
out["table"] = rows

# ---------------------------------------------------------------- fig: element fields
th = np.linspace(0, np.pi, 181)
ph0 = np.zeros_like(th)
d, hx, rg = ideal.dipole(), ideal.helix(), ideal.ring(R=R_RET)
tot = Elements.concat([d, hx, rg])
ff = {n: far_field(e, th, ph0, k) for n, e in [("dipole", d), ("coil", hx), ("ret", rg), ("total", tot)]}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
a = ax[0]
for n, lab in [("dipole", "dipole"), ("coil", "toroidal coil"), ("total", "sum of all parts")]:
    a.plot(np.degrees(th), -ff[n][0].real / H_D, color=COL[n], label=lab, lw=2.4 if n == "total" else 2)
a.set_title(r"$\theta$-polarized field (in phase with $I$)")
a.set_xlabel(r"$\theta$  (deg, from the dipole axis)"); a.set_ylabel(r"signed far field,  $-$Re $N_\theta$ / ($I\,h$)")
a.set_xticks(range(0, 181, 30)); a.legend(loc="upper right", fontsize=8.5); a.axhline(0, color=AXIS, lw=1)
a.text(90, -0.55, "the coil field is equal\nand opposite to the dipole's", color=INK2, fontsize=9, ha="center")
a = ax[1]
for n, lab in [("coil", r"coil: one turn around the ring ($m_z>0$)"), ("ret", r"return wire ($m_z<0$)"),
               ("total", "sum of all parts")]:
    a.plot(np.degrees(th), ff[n][1].imag / H_D, color=COL[n], label=lab, lw=2.4 if n == "total" else 2)
a.set_title(r"$\varphi$-polarized field (in quadrature with $I$)")
a.set_xlabel(r"$\theta$  (deg)"); a.set_ylabel(r"Im $N_\varphi$ / ($I\,h$)")
a.set_xticks(range(0, 181, 30)); a.legend(loc="upper right", fontsize=8.5); a.axhline(0, color=AXIS, lw=1); a.set_ylim(-6, 9)
save(fig, "m0_element_fields")
out["equator"] = {n: dict(Nth=[ff[n][0][90].real / H_D, ff[n][0][90].imag / H_D],
                          Nph=[ff[n][1][90].real / H_D, ff[n][1][90].imag / H_D]) for n in ff}

# ---------------------------------------------------------------- fig: 3D patterns
T, P = np.meshgrid(np.linspace(0, np.pi, 61), np.linspace(0, 2 * np.pi, 121), indexing="ij")
pat = [("dipole alone", ideal.dipole(), COL["dipole"]),
       ("toroidal coil alone", ideal.loops(), COL["coil"]),
       ("dipole + coil (anapole)", ideal.anapole_loops(), COL["total"]),
       ("dipole + reversed coil", ideal.anapole_loops(sense=-1), COL["rev"])]
fig = plt.figure(figsize=(13, 3.9))
for i, (name, e, c) in enumerate(pat):
    Nt, Np = far_field(e, T, P, k)
    A = np.sqrt(abs(Nt) ** 2 + abs(Np) ** 2) / H_D
    a = fig.add_subplot(1, 4, i + 1, projection="3d")
    X, Y, Z = A * np.sin(T) * np.cos(P), A * np.sin(T) * np.sin(P), A * np.cos(T)
    a.plot_surface(X, Y, Z, color=c, alpha=0.85, linewidth=0.15, edgecolor=SURF, rstride=2, cstride=4, shade=True)
    a.plot([0, 0], [0, 0], [-1.6, 1.6], color=MUTED, lw=0.8)
    a.set_xlim(-2, 2); a.set_ylim(-2, 2); a.set_zlim(-1.6, 1.6); a.set_box_aspect((1, 1, 0.8)); a.set_axis_off()
    a.view_init(22, -55)
    a.set_title("%s\npeak field %.3f × dipole (%+.1f dB)" % (name, A.max(), 20 * np.log10(A.max())),
                fontsize=8.5, loc="center", fontweight="normal")
fig.subplots_adjust(wspace=-0.05)
save(fig, "m0_patterns_3d")

# ---------------------------------------------------------------- fig: residual pattern
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
Tg, Pg = np.meshgrid(np.linspace(0, np.pi, 91), np.linspace(0, 2 * np.pi, 181), indexing="ij")
for a, (h, lab) in zip(ax, [(H_D, "h = 29.80 mm (design)"), (0.02946, "h = 29.46 mm (tuned)")]):
    Nt, Np = far_field(ideal.anapole_loops(h=h), Tg, Pg, k)
    D = 20 * np.log10(np.sqrt(abs(Nt) ** 2 + abs(Np) ** 2) / H_D + 1e-12)
    im = a.pcolormesh(np.degrees(Pg), np.degrees(Tg), D, cmap="Blues", vmin=-80, vmax=-30, shading="auto")
    a.set_title("%s — max %.1f dB" % (lab, D.max())); a.set_xlabel(r"$\varphi$ (deg)"); a.grid(False)
    a.set_xticks(range(0, 361, 60))
ax[0].set_ylabel(r"$\theta$ (deg)"); ax[0].invert_yaxis()
cb = fig.colorbar(im, ax=ax, pad=0.02); cb.set_label("residual far field (dB, bare dipole = 0)")
save(fig, "m0_residual_map")

# ---------------------------------------------------------------- sweeps
hs = np.linspace(0.026, 0.034, 161)
coil = ideal.loops()
T, P, W = sphere_grid(41, 72)
Nc = far_field(coil, T, P, k); Nd1 = far_field(ideal.dipole(h=1.0, n=1), T, P, k)   # dipole ~ linear in h (kh<<1)
res_h = []
for h in hs:
    Nd = far_field(ideal.dipole(h=h), T, P, k)
    U = abs(Nd[0] + Nc[0]) ** 2 + abs(Nd[1] + Nc[1]) ** 2
    Ud = abs(far_field(ideal.dipole(), T, P, k)[0]) ** 2
    res_h.append(10 * np.log10((U * W).sum() / (Ud * W).sum()))
res_h = np.array(res_h)
fs = np.linspace(90e6, 120e6, 301)
res_f = []
an = ideal.anapole_loops(); dd = ideal.dipole()
for f in fs:
    kk = k_of(f)
    res_f.append(10 * np.log10(radiated_power(an, kk, 41, 72) / radiated_power(dd, kk, 41, 72)))
res_f = np.array(res_f)
out["h_null_mm"] = float(hs[res_h.argmin()] * 1e3); out["f_null_MHz"] = float(fs[res_f.argmin()] / 1e6)
out["res_at_f0_db"] = float(np.interp(F0, fs, res_f))
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
ax[0].plot(hs * 1e3, res_h, color=COL["total"], label="exact integration")
ax[0].plot(hs * 1e3, 20 * np.log10(abs(1 - H_T / hs * 1) * hs / H_D + 1e-9), color=COL["theory"], ls="--",
           label=r"2nd-order theory $(h-h_T)^2$")
ax[0].axvline(H_D * 1e3, color=AXIS, lw=1); ax[0].text(H_D * 1e3 + .05, -8, "design h", color=INK2, fontsize=9)
ax[0].set_xlabel("dipole length h (mm)"); ax[0].set_ylabel("radiated power / bare 29.8 mm dipole (dB)")
ax[0].set_title("Null vs dipole length (f = 104.2 MHz)"); ax[0].legend(loc="lower left")
ax[1].plot(fs / 1e6, res_f, color=COL["total"], label="exact integration")
ax[1].plot(fs / 1e6, 20 * np.log10(abs(1 - (fs / F0) ** 2) + 1e-9), color=COL["theory"], ls="--",
           label=r"2nd-order theory $[1-(f/f_0)^2]^2$")
ax[1].axvline(F0 / 1e6, color=AXIS, lw=1); ax[1].text(F0 / 1e6 + .2, -8, "$f_0$", color=INK2, fontsize=9)
ax[1].set_xlabel("frequency (MHz)"); ax[1].set_title("Null vs frequency (h = 29.8 mm)"); ax[1].legend(loc="lower left")
ax[0].set_ylim(-75, 0)
save(fig, "m0_sweeps")

# ---------------------------------------------------------------- near field maps
dN, cN = ideal.dipole(n=20), ideal.loops(n=72)
aN = Elements.concat([dN, cN])
ext = 0.6; n = 171
x = np.linspace(-ext, ext, n); X, Zz = np.meshgrid(x, x, indexing="xy")
planes = {"xz": np.stack([X, 0 * X + 1e-4, Zz], -1), "xy": np.stack([X, Zz, 0 * X + 1e-4], -1)}
nf = {}
for pl, pts in planes.items():
    for nme, e in [("dipole", dN), ("coil", cN), ("total", aN)]:
        E, H = near_field(e, pts.reshape(-1, 3), k)
        nf[pl, nme] = (np.linalg.norm(E, axis=1).reshape(n, n), np.linalg.norm(H, axis=1).reshape(n, n),
                       E.reshape(n, n, 3), H.reshape(n, n, 3))
tt = np.linspace(0, 2 * np.pi, 100)


def outline(a, pl):
    if pl == "xz":
        for s in (-1, 1):
            a.plot(s * R_T + A_L * np.cos(tt), A_L * np.sin(tt), color="white", lw=0.8)
        a.plot([0, 0], [-H_D / 2, H_D / 2], color="white", lw=1.5)
    else:
        for i in range(N_L):
            p = 2 * np.pi * i / N_L
            a.plot([(R_T - A_L) * np.cos(p), (R_T + A_L) * np.cos(p)], [(R_T - A_L) * np.sin(p), (R_T + A_L) * np.sin(p)],
                   color="white", lw=0.8)
    a.plot((R_T + A_L) * np.cos(tt), (R_T + A_L) * np.sin(tt), color="white", lw=0.6, ls=":")


for comp, idx, unit, vmin, vmax, pl, ttl in [("E", 0, "V/m", 1e-1, 1e5, "xz", "Electric field |E|, meridian plane (x–z)"),
                                             ("H", 1, "A/m", 1e-4, 1e2, "xy", "Magnetic field |H|, equatorial plane (x–y)")]:
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
    for a, (nme, lab) in zip(ax, [("dipole", "dipole alone"), ("coil", "toroidal coil alone"), ("total", "dipole + coil (anapole)")]):
        im = a.pcolormesh(x, x, np.clip(nf[pl, nme][idx], vmin, vmax), norm=LogNorm(vmin, vmax), cmap="magma", shading="auto")
        outline(a, pl); a.set_aspect("equal"); a.grid(False); a.set_title(lab); a.set_xlabel("x (m)")
    ax[0].set_ylabel("z (m)" if pl == "xz" else "y (m)")
    cb = fig.colorbar(im, ax=ax, pad=0.015, shrink=0.9); cb.set_label("|%s| (%s) for I = 1 A" % (comp, unit))
    fig.suptitle(ttl + " — dotted circle: sphere r = R + a that encloses the structure", x=0.07, ha="left", fontsize=11, color=INK)
    save(fig, "m0_near_%s_%s" % (comp, pl))

# ---------------------------------------------------------------- radial decay + power flux
rr = np.geomspace(0.02, 30, 140)
dirs = {"equator (θ=90°, φ=15°)": (np.pi / 2, np.radians(15)), "θ=45°, φ=15°": (np.pi / 4, np.radians(15))}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.1))
a = ax[0]
thd, phd = dirs["equator (θ=90°, φ=15°)"]
u = np.array([np.sin(thd) * np.cos(phd), np.sin(thd) * np.sin(phd), np.cos(thd)])
dec = {}
for nme, e, lab in [("dipole", dN, "dipole alone"), ("coil", cN, "coil alone"), ("total", aN, "dipole + coil (anapole)")]:
    E, H = near_field(e, rr[:, None] * u[None], k)
    dec[nme] = np.linalg.norm(E, axis=1)
    a.loglog(rr, dec[nme], color=COL[nme], label=lab)
a.axvline(R_T + A_L, color=AXIS, lw=1); a.text(R_T + A_L * 1.15, 3e4, "r = R + a", color=INK2, fontsize=9)
a.axvline(LAM0 / (2 * np.pi), color=AXIS, lw=1, ls=":"); a.text(LAM0 / (2 * np.pi) * 1.08, 3e3, r"r = $\lambda/2\pi$", color=INK2, fontsize=9)
a.set_xlabel("distance from the centre r (m)"); a.set_ylabel("|E| (V/m) for I = 1 A")
a.set_title("Field along the equator"); a.legend(loc="lower left")
a.set_ylim(1e-5, 1e6)
out["decay"] = dict(r=rr.tolist(), **{n: v.tolist() for n, v in dec.items()})
# power through spheres
a = ax[1]
T, P, W = sphere_grid(40, 96)
nrm = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], -1).reshape(-1, 3)
rs = np.geomspace(0.25, 20, 16)
flux = {}
for nme, e, lab in [("dipole", dN, "dipole alone"), ("coil", cN, "coil alone"), ("total", aN, "dipole + coil (anapole)")]:
    Pr, Qr = [], []
    for r0 in rs:
        E, H = near_field(e, r0 * nrm, k)
        S = 0.5 * np.cross(E, np.conj(H))
        Sn = (S * nrm).sum(-1) * W.ravel() * r0 ** 2
        Pr.append(Sn.real.sum()); Qr.append(abs(Sn.imag.sum()))
    flux[nme] = (np.array(Pr), np.array(Qr))
    a.loglog(rs, np.abs(flux[nme][0]), color=COL[nme], label=lab + " — active power")
    a.loglog(rs, flux[nme][1], color=COL[nme], ls=":", lw=1.4)
a.plot([], [], color=MUTED, ls=":", lw=1.4, label="reactive power (same colours)")
a.set_xlabel("radius of the sphere r (m)"); a.set_ylabel("power through the sphere (W) for I = 1 A")
a.set_title("Power crossing a sphere of radius r"); a.legend(loc="lower left", fontsize=8)
save(fig, "m0_decay_flux")
out["flux"] = dict(r=rs.tolist(), **{n: dict(P=v[0].tolist(), Q=v[1].tolist()) for n, v in flux.items()})
out["suppression_outside_db"] = float(20 * np.log10(dec["total"][-1] / dec["dipole"][-1]))
i1 = np.argmin(abs(rr - 0.3))
out["suppression_at_0.3m_db"] = float(20 * np.log10(dec["total"][i1] / dec["dipole"][i1]))

# stored energy inside / outside (reactive) not computed here
json.dump(out, open("out/model0.json", "w"), indent=1)
print(json.dumps({k_: v for k_, v in out.items() if k_ not in ("table", "decay", "flux")}, indent=1))
