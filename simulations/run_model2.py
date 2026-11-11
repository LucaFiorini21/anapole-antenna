"""Model 2 - NEC-2, single generator: dipole, toroidal winding and return
ring all in series (feeding option 1 of the report), as one continuous wire."""
import json
import numpy as np
from anapole.core import *
from anapole.models import *
from anapole.analysis import *
from anapole.style import *

k = K0
out = {}
Pd, p0d, Zd, rd = bare_dipole()


def run(f=F0, **kw):
    st = model_series(**kw)
    r = st.solve(f)
    If = r.current(*st.feed)
    Z = 1 / If
    r.I = r.I / If                      # unit current at the generator
    return st, r, Z


# ---- design frequency, several segmentations ------------------------------
rows = []
for kw in [dict(n_t=12), dict(n_t=16), dict(n_t=24), dict(n_t=32), dict(n_t=16, dz=0.004), dict(n_t=16, rad=0.0005),
           dict(n_t=16, sense=-1), dict(n_t=16, with_return=False)]:
    st, r, Z = run(**kw)
    el = r.elements(); F = ff(el, k)
    part = np.array(r.part)
    Ic = r.I[part == "coil"]
    rows.append(dict(kw=kw, nseg=len(r.I), Z=Z, P_db=db(prad_ff(F, k) / Pd), wire_m=st.wire_length("coil") + st.wire_length("return"),
                     coil_Imin=abs(Ic).min(), coil_Imax=abs(Ic).max(), I_upper_arm=abs(r.I[part == "dipole"][-1]),
                     parts=part_table(el, k)))
    print("%-34s nseg=%3d Z=%7.1f%+8.0fj  P/Pdip=%+6.1f dB  coil |I| %.2f..%.2f  upper-arm |I|=%.2f" %
          (kw, len(r.I), Z.real, Z.imag, rows[-1]["P_db"], abs(Ic).min(), abs(Ic).max(), rows[-1]["I_upper_arm"]))
out["f0"] = rows

# ---- current along the wire ----------------------------------------------
st, r, Z = run(n_t=24)
part = np.array(r.part)
order = np.where(np.isin(part, ["dipole", "leads", "coil", "return"]))[0]
# wire order in the model: lower arm, upper arm, plates..., L1, helix, J1, ring, L2
tags = r.tag
seq = [st.dip_ref[0], st.feed[0]] + st.helix_tags + [st.helix_tags[-1] + 1] + st.ret_tags + [st.ret_tags[-1] + 1, st.dip_ref[0] + 1]
idx = np.concatenate([np.where(tags == t)[0] for t in seq])
s = np.cumsum(np.linalg.norm(r.L[idx], axis=1)); s -= s[0]
I = r.I[idx]; pp = part[idx]
fig, ax = plt.subplots(2, 1, figsize=(11, 5.2), sharex=True, gridspec_kw=dict(hspace=0.12))
names = {"dipole": ("dipole", COL["dipole"]), "leads": ("leads", COL["leads"]), "coil": ("toroidal winding", COL["coil"]), "return": ("return ring", COL["ret"])}
for pn, (lab, c) in names.items():
    m = pp == pn
    ax[0].plot(np.where(m, s, np.nan) / LAM0, np.where(m, abs(I), np.nan), color=c, label=lab)
    ax[1].plot(np.where(m, s, np.nan) / LAM0, np.where(m, np.degrees(np.angle(I)), np.nan), color=c)
ax[0].axhline(1, color=MUTED, ls="--", lw=1); ax[0].text(0.52, 2.2, "dashed: what the theory assumes\n(|I| = 1, phase 0° everywhere)", color=INK2, fontsize=9)
ax[0].set_ylabel("|I| / I at the generator"); ax[0].legend(loc="upper right", ncol=4, fontsize=8.5)
ax[0].set_title("Series antenna, one generator: current along the %.1f m of wire at 104.2 MHz" % s[-1])
ax[1].set_ylabel("phase of I (deg)"); ax[1].set_yticks([-180, -90, 0, 90, 180]); ax[1].axhline(0, color=MUTED, ls="--", lw=1)
ax[1].set_xlabel("position along the wire, from the bottom plate to the top plate (wavelengths)")
save(fig, "m2_current_along_wire")
out["wire_len_m"] = float(s[-1]); out["wire_len_lambda"] = float(s[-1] / LAM0)

# ---- frequency sweep -------------------------------------------------------
fs = np.geomspace(3e6, 130e6, 150)
sw = []
for f in fs:
    kk = k_of(f)
    st, r, Z = run(f, n_t=16)
    el = r.elements()
    part = np.array(r.part)
    Iw = r.I[np.isin(part, ["coil", "return"])]
    # bare dipole with plates at this frequency (same generator current)
    Pdf = bare_dipole(f)[0]
    P = radiated_power(el, kk, 31, 60)
    sw.append((db(P / Pdf), np.std(Iw) / abs(np.mean(Iw)), abs(Iw).min(), abs(Iw).max(), Z.real, Z.imag))
sw = np.array(sw)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
a = ax[0]
a.semilogx(fs / 1e6, sw[:, 0], color=COL["total"], label="NEC-2, series antenna")
a.semilogx(fs / 1e6, 20 * np.log10(abs(1 - (fs / F0) ** 2) + 1e-6), color=COL["theory"], ls="--", label=r"theory, uniform current: $[1-(f/f_0)^2]^2$")
a.axvline(F0 / 1e6, color=AXIS, lw=1); a.text(F0 / 1e6 * 0.97, 32, "$f_0$", ha="right", color=INK2)
a.set_ylim(-60, 40); a.set_xlabel("frequency (MHz)"); a.set_ylabel("radiated power vs bare dipole, same current (dB)")
a.set_title("The predicted null never appears"); a.legend(loc="lower left", fontsize=8.5)
a = ax[1]
a.loglog(fs / 1e6, np.minimum(sw[:, 3] / np.maximum(sw[:, 2], 1e-3), 1e3), color=COL["coil"])
a.axvline(F0 / 1e6, color=AXIS, lw=1); a.axhline(1, color=MUTED, ls="--", lw=1)
a.set_xlabel("frequency (MHz)"); a.set_ylabel("max |I| / min |I| along winding + ring")
a.set_title("The current is never uniform along the wire"); a.set_ylim(0.8, 1.5e3)
a.text(3.2, 1.12, "uniform current = 1", color=INK2, fontsize=9)
save(fig, "m2_frequency_sweep")
out["sweep"] = dict(f=fs, P_db=sw[:, 0], nonuni=sw[:, 1])
i10 = np.where(sw[:, 3] / np.maximum(sw[:, 2], 1e-3) > 1.1)[0]
out["f_uniform_within_10pct_MHz"] = float(fs[i10[0]] / 1e6) if len(i10) else None
print("current non-uniform by >10%% above %.1f MHz" % out["f_uniform_within_10pct_MHz"])

st = model_series(n_t=16)
open("nec/model2_series_single_generator.nec", "w").write(st.deck(F0, "Anapole attempt - dipole + toroidal winding + return ring in series, one generator"))
json.dump(jsonable(out), open("out/model2.json", "w"), indent=1)
