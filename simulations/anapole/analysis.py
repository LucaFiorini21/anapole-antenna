"""Analysis helpers shared by the NEC-2 runs."""
import numpy as np
from .core import *
from .models import model_loops, solve_ports

_T, _P, _W = sphere_grid(41, 80)


def ff(el, k):
    """far field on the sphere grid, stacked [N_theta, N_phi] (2, nth, nph)."""
    return np.stack(far_field(el, _T, _P, k))


def prad_ff(F, k):
    return float(ETA * k ** 2 / (32 * np.pi ** 2) * ((abs(F) ** 2).sum(0) * _W).sum())


def peak_ff(F):
    return float(np.sqrt((abs(F) ** 2).sum(0)).max())


def port_basis(st, f, lossy=False):
    """For every port group g: Result with unit current in g's reference port
    and zero current in the reference ports of the other groups."""
    groups = list(st.ports)
    res = {}
    for g in groups:
        r, V, Y = solve_ports(st, f, {x: (1.0 if x == g else 0.0) for x in groups}, lossy)
        res[g] = (r, V)
    return res


def best_combination(Fs, fixed=0):
    """complex weights w (w[fixed]=1) minimising the radiated power of
    sum_j w_j F_j  (least squares on the sphere)."""
    sw = np.sqrt(_W)[None]
    cols = [(F * sw).ravel() for F in Fs]
    A = np.stack([c for j, c in enumerate(cols) if j != fixed], 1)
    b = -cols[fixed]
    x = np.linalg.lstsq(A, b, rcond=None)[0]
    w = np.insert(x, fixed, 1.0)
    return w


_ref = {}


def bare_dipole(f=F0, **kw):
    """NEC reference: the dipole with end plates alone, unit feed current.
    returns (P_rad, p0z, Z_in, Result)."""
    key = (f, tuple(sorted(kw.items())))
    if key not in _ref:
        st = model_loops(with_loops=False, **kw)
        r, V, Y = solve_ports(st, f, {"dipole": 1.0})
        el = r.elements()
        _ref[key] = (prad_ff(ff(el, r.k), r.k), moments(el, r.k)["p0"][2], V["dipole"], r)
    return _ref[key]


def db(x):
    return 10 * np.log10(np.maximum(x, 1e-30))


def part_table(el, k):
    """moments of every part: p0z (mm), m_z (cm^2), k^2 T_z (mm), complex."""
    rows = {}
    for n in dict.fromkeys(el.names.values()):
        e = el.sel(n)
        if len(e.I) == 0:
            continue
        m = moments(e, k)
        rows[n] = dict(p0z_mm=complex(m["p0"][2] * 1e3), mz_cm2=complex(m["m"][2] * 1e4),
                       k2Tz_mm=complex(k * k * m["T"][2] * 1e3),
                       Imin=float(abs(e.I).min()), Imax=float(abs(e.I).max()))
    return rows


def jsonable(o):
    if isinstance(o, dict):
        return {str(k_): jsonable(v) for k_, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (complex, np.complexfloating)):
        return [float(np.real(o)), float(np.imag(o))]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    return o
