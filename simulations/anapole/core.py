"""Core tools for the anapole-antenna study.

* constants and design values (from the repository README)
* Structure: a list of straight wires, grouped in named parts, that can be
  - solved with NEC-2 (PyNEC / necpp, method of moments), or
  - exported as a .nec deck for 4nec2 / xnec2c
* fields(): near and far field of an arbitrary set of current elements,
  from the exact Hertzian-dipole expressions (no multipole expansion).

Time convention exp(+j w t), SI units.
"""
import numpy as np

C0 = 299792458.0
MU0 = 4e-7 * np.pi
ETA = MU0 * C0
EPS0 = 1.0 / (MU0 * C0 ** 2)

# ---- design values (README, sections 1-2) -------------------------------
F0 = 104.2e6          # Hz
R_T = 0.15            # ring radius R  [m]
A_L = 0.047           # loop radius a  [m]
N_L = 12              # number of loops
H_D = 0.0298          # dipole length h [m]
WIRE_R = 0.001        # wire radius (2 mm copper wire) [m]
SIGMA_CU = 5.8e7


def k_of(f):
    return 2 * np.pi * f / C0


K0 = k_of(F0)
LAM0 = C0 / F0
H_T = K0 ** 2 * N_L * np.pi * A_L ** 2 * R_T / 2   # toroidal effective length


# =========================================================================
#  Current elements and fields
# =========================================================================
class Elements:
    """A set of short straight current filaments.
    r  : (K,3) centres,  L : (K,3) vector lengths,  I : (K,) complex currents
    part: (K,) integer label of the part each element belongs to."""

    def __init__(self, r, L, I, part=None, names=None):
        self.r = np.asarray(r, float)
        self.L = np.asarray(L, float)
        self.I = np.asarray(I, complex)
        self.part = np.zeros(len(self.r), int) if part is None else np.asarray(part, int)
        self.names = names or {0: "all"}

    def sel(self, *names):
        ids = [i for i, n in self.names.items() if n in names]
        m = np.isin(self.part, ids)
        return Elements(self.r[m], self.L[m], self.I[m], self.part[m], self.names)

    def scaled(self, s):
        return Elements(self.r, self.L, self.I * s, self.part, self.names)

    @staticmethod
    def concat(lst):
        names = {}
        r, L, I, p = [], [], [], []
        off = 0
        for e in lst:
            mp = {}
            for i, n in e.names.items():
                if n in names.values():
                    mp[i] = [k for k, v in names.items() if v == n][0]
                else:
                    mp[i] = off
                    names[off] = n
                    off += 1
            r.append(e.r); L.append(e.L); I.append(e.I)
            p.append(np.array([mp[i] for i in e.part], int))
        return Elements(np.vstack(r), np.vstack(L), np.concatenate(I), np.concatenate(p), names)

    @property
    def p(self):  # current moments I*L  (K,3)
        return self.I[:, None] * self.L


def sph_dirs(theta, phi):
    """unit vectors n, theta_hat, phi_hat for arrays theta, phi (radians)."""
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    n = np.stack([st * cp, st * sp, ct], -1)
    th = np.stack([ct * cp, ct * sp, -st], -1)
    ph = np.stack([-sp, cp, 0 * st], -1)
    return n, th, ph


def radiation_vector(el, n, k):
    """N(n) = sum I L exp(+jk n.r')  ->  (...,3) complex."""
    n2 = n.reshape(-1, 3)
    ph = np.exp(1j * k * n2 @ el.r.T)               # (M,K)
    # finite-length correction of each straight element (sinc), exact for uniform I
    arg = 0.5 * k * n2 @ el.L.T
    ph = ph * np.sinc(arg / np.pi)
    N = ph @ el.p                                    # (M,3)
    return N.reshape(n.shape)


def far_field(el, theta, phi, k):
    """returns N_theta, N_phi (transverse radiation vector, in A*m).
    E_far * r * exp(jkr) = -j k eta/(4 pi) * N_perp."""
    n, th, ph = sph_dirs(theta, phi)
    N = radiation_vector(el, n, k)
    return (N * th).sum(-1), (N * ph).sum(-1)


def sphere_grid(nth=91, nph=180):
    """Gauss-Legendre in cos(theta) x uniform in phi; weights sum to 4 pi."""
    x, w = np.polynomial.legendre.leggauss(nth)
    theta = np.arccos(x)
    phi = np.arange(nph) * 2 * np.pi / nph
    T, P = np.meshgrid(theta, phi, indexing="ij")
    W = np.outer(w, np.full(nph, 2 * np.pi / nph))
    return T, P, W


def radiated_power(el, k, nth=61, nph=120):
    """P_rad [W] for the given (peak-phasor) currents, by integration of the
    far field over the sphere:  U = eta k^2/(32 pi^2) |N_perp|^2."""
    T, P, W = sphere_grid(nth, nph)
    Nt, Np = far_field(el, T, P, k)
    U = ETA * k ** 2 / (32 * np.pi ** 2) * (abs(Nt) ** 2 + abs(Np) ** 2)
    return float((U * W).sum())


def near_field(el, pts, k, chunk=4000):
    """Exact E, H at points pts (M,3) from the Hertzian-dipole fields of all
    elements (includes the charges implied by current discontinuities)."""
    pts = np.asarray(pts, float).reshape(-1, 3)
    E = np.zeros((len(pts), 3), complex)
    H = np.zeros((len(pts), 3), complex)
    p = el.p
    for s in range(0, len(pts), chunk):
        d = pts[s:s + chunk, None, :] - el.r[None, :, :]        # (m,K,3)
        r = np.linalg.norm(d, axis=-1)
        rh = d / r[..., None]
        kr = k * r
        g = np.exp(-1j * kr) / (4 * np.pi * r)
        pr = (p[None] * rh).sum(-1)                               # p . r_hat
        pperp = p[None] - pr[..., None] * rh
        a = 1 + 1 / (1j * kr)
        Et = -1j * k * ETA * (g * (a - 1 / kr ** 2))[..., None] * pperp
        Er = (2 * ETA * g / r * a * pr)[..., None] * rh
        E[s:s + chunk] = (Et + Er).sum(1)
        H[s:s + chunk] = (1j * k * (g * a)[..., None] * np.cross(p[None], rh)).sum(1)
    return E, H


def moments(el, k):
    """Cartesian multipole moments of the current set (origin = centre).
    Returns dict: p0 = int J dV (A m), m (A m^2), T (A m^3)."""
    r, p = el.r, el.p
    p0 = p.sum(0)
    m = 0.5 * np.cross(r, p).sum(0)
    rp = (r * p).sum(-1)
    r2 = (r * r).sum(-1)
    T = 0.1 * (rp[:, None] * r - 2 * r2[:, None] * p).sum(0)
    return dict(p0=p0, m=m, T=T)


# =========================================================================
#  Wire structures for NEC-2
# =========================================================================
class Structure:
    def __init__(self):
        self.wires = []      # dict(tag,nseg,p1,p2,rad,part)
        self.sources = []    # (tag, seg(1-based), complex V)
        self.loads = []      # (tag, seg, R, L, C) series lumped loads
        self._tag = 0

    def wire(self, p1, p2, nseg, part, rad=WIRE_R):
        self._tag += 1
        self.wires.append(dict(tag=self._tag, nseg=int(nseg), p1=np.array(p1, float),
                               p2=np.array(p2, float), rad=rad, part=part))
        return self._tag

    def polyline(self, pts, part, rad=WIRE_R, seg_per=1):
        """one wire per polyline edge; returns list of tags."""
        return [self.wire(pts[i], pts[i + 1], seg_per, part, rad) for i in range(len(pts) - 1)]

    def load(self, tag, seg, R=0.0, L=0.0, C=0.0):
        self.loads.append((tag, seg, R, L, C))

    def source(self, tag, seg, V=1.0):
        self.sources.append((tag, seg, complex(V)))

    # ---- segment table --------------------------------------------------
    def segments(self):
        r, L, part, tag, seg = [], [], [], [], []
        for w in self.wires:
            d = (w["p2"] - w["p1"]) / w["nseg"]
            for i in range(w["nseg"]):
                r.append(w["p1"] + (i + 0.5) * d); L.append(d)
                part.append(w["part"]); tag.append(w["tag"]); seg.append(i + 1)
        return np.array(r), np.array(L), part, np.array(tag), np.array(seg)

    def seg_index(self, tag, seg):
        _, _, _, t, s = self.segments()
        return int(np.where((t == tag) & (s == seg))[0][0])

    def wire_length(self, part=None):
        return sum(np.linalg.norm(w["p2"] - w["p1"]) for w in self.wires
                   if part is None or w["part"] == part)

    # ---- NEC deck -------------------------------------------------------
    def deck(self, f_hz, title="anapole antenna", lossy=False, rp=True):
        s = [f"CM {title}", "CM generated by anapole/core.py (units: m, MHz)", "CE"]
        for w in self.wires:
            a, b = w["p1"], w["p2"]
            s.append("GW %d %d %.7f %.7f %.7f %.7f %.7f %.7f %.5f" %
                     (w["tag"], w["nseg"], *a, *b, w["rad"]))
        s += ["GE 0", "EK"]
        if lossy:
            s.append("LD 5 0 0 0 %.4e" % SIGMA_CU)
        for t, sg, R, L, C in self.loads:
            s.append("LD 0 %d %d %d %.6e %.6e %.6e" % (t, sg, sg, R, L, C))
        s.append("FR 0 1 0 0 %.6f 0" % (f_hz / 1e6))
        for t, sg, V in self.sources:
            s.append("EX 0 %d %d 0 %.8e %.8e" % (t, sg, V.real, V.imag))
        if rp:
            s.append("RP 0 37 73 1000 0 0 5 5")
        s.append("EN")
        return "\n".join(s) + "\n"

    # ---- NEC solve ------------------------------------------------------
    def solve(self, f_hz, lossy=False, sources=None, ext_kernel=True):
        """Run NEC-2. Returns Result with segment currents (peak phasors)."""
        from PyNEC import nec_context
        ctx = nec_context()
        g = ctx.get_geometry()
        for w in self.wires:
            a, b = w["p1"], w["p2"]
            g.wire(w["tag"], w["nseg"], a[0], a[1], a[2], b[0], b[1], b[2], w["rad"], 1.0, 1.0)
        ctx.geometry_complete(0)
        ctx.set_extended_thin_wire_kernel(bool(ext_kernel))
        if lossy:
            ctx.ld_card(5, 0, 0, 0, SIGMA_CU, 0.0, 0.0)
        for t, sg, R, L, C in self.loads:
            ctx.ld_card(0, t, sg, sg, R, L, C)
        ctx.fr_card(0, 1, f_hz / 1e6, 0.0)
        src = self.sources if sources is None else sources
        for t, sg, V in src:
            ctx.ex_card(0, t, sg, 0, V.real, V.imag, 0, 0, 0, 0)
        ctx.xq_card(0)
        sc = ctx.get_structure_currents(0)
        I = np.array(sc.get_current())
        ip = ctx.get_input_parameters(0)
        Vs = np.array(ip.get_voltage()); Is = np.array(ip.get_current())
        return Result(self, f_hz, I, Vs, Is, list(src))


class Result:
    def __init__(self, st, f, I, Vs, Is, src):
        self.st, self.f, self.k = st, f, k_of(f)
        self.I, self.Vs, self.Is, self.src = I, Vs, Is, src
        self.r, self.L, self.part, self.tag, self.seg = st.segments()

    @property
    def p_in(self):
        return float(0.5 * np.real(self.Vs * np.conj(self.Is)).sum())

    def current(self, tag, seg):
        return self.I[self.st.seg_index(tag, seg)]

    def elements(self):
        names = {}
        ids = []
        for p in self.part:
            if p not in names.values():
                names[len(names)] = p
            ids.append([k for k, v in names.items() if v == p][0])
        return Elements(self.r, self.L, self.I, ids, names)

    def p_loss(self):
        """ohmic loss with skin-effect surface resistance of copper."""
        rs = np.sqrt(np.pi * self.f * MU0 / SIGMA_CU)
        rad = np.concatenate([[w["rad"]] * w["nseg"] for w in self.st.wires])
        l = np.linalg.norm(self.L, axis=1)
        return float(0.5 * (abs(self.I) ** 2 * rs * l / (2 * np.pi * rad)).sum())
