"""Model 0 - prescribed (ideal) currents: the assumptions of the report.
Uniform current I in the dipole (point charges at its ends) and the same
current I in every turn of the coil."""
import numpy as np
from .core import Elements, R_T, A_L, N_L, H_D


def dipole(h=H_D, I=1.0, n=40):
    z = (np.arange(n) + 0.5) / n * h - h / 2
    r = np.stack([0 * z, 0 * z, z], -1)
    L = np.tile([0, 0, h / n], (n, 1))
    return Elements(r, L, np.full(n, I, complex), names={0: "dipole"})


def _curve(pts, I, name):
    """elements from a polyline pts (n+1,3) (centre + vector length)."""
    r = 0.5 * (pts[1:] + pts[:-1]); L = pts[1:] - pts[:-1]
    return Elements(r, L, np.full(len(r), I, complex), names={0: name})


def loops(R=R_T, a=A_L, N=N_L, I=1.0, n=240, sense=+1):
    """N separate circular loops in meridian planes.
    sense=+1: current along +z on the OUTER side (m_i = -I A phi_hat)."""
    psi = np.linspace(0, 2 * np.pi, n + 1) * sense
    out = []
    for i in range(N):
        ph = 2 * np.pi * i / N
        rho = R + a * np.cos(psi)
        pts = np.stack([rho * np.cos(ph), rho * np.sin(ph), a * np.sin(psi)], -1)
        out.append(_curve(pts, I, "coil"))
    return Elements.concat(out)


def helix(R=R_T, a=A_L, N=N_L, I=1.0, n=240, sense=+1, adv=+1, psi0=0.0):
    """continuous toroidal winding: N turns while advancing once in phi."""
    t = np.linspace(0, 1, N * n + 1)
    psi = psi0 + 2 * np.pi * N * t * sense
    ph = 2 * np.pi * t * adv
    rho = R + a * np.cos(psi)
    pts = np.stack([rho * np.cos(ph), rho * np.sin(ph), a * np.sin(psi)], -1)
    return _curve(pts, I, "coil")


def ring(R=R_T, I=1.0, n=720, adv=-1, z=0.0):
    ph = np.linspace(0, 2 * np.pi, n + 1) * adv
    pts = np.stack([R * np.cos(ph), R * np.sin(ph), 0 * ph + z], -1)
    return _curve(pts, I, "return")


def anapole_loops(h=H_D, sense=+1, **kw):
    return Elements.concat([dipole(h), loops(sense=sense, **kw)])


def anapole_helix(h=H_D, sense=+1, with_return=True, **kw):
    parts = [dipole(h), helix(sense=sense, **kw)]
    if with_return:
        parts.append(ring())
    return Elements.concat(parts)
