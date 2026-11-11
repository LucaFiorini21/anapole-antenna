"""Wire models for NEC-2.

model_loops : dipole with end plates + N independent loops (each with its own
              feed point(s)).  Tests the theory with real current distributions.
model_series: dipole + continuous toroidal winding + return ring, all in
              series, one generator.  The 'as built' single-feed antenna.
"""
import numpy as np
from .core import Structure, R_T, A_L, N_L, H_D, WIRE_R


def _eq_area(n):
    """radius scale so that a regular n-gon has the area of the circle."""
    return np.sqrt(2 * np.pi / (n * np.sin(2 * np.pi / n)))


def add_dipole(st, h, nd=5, r_p=0.03, n_sp=6, sp_seg=3, zgap=None, rad=WIRE_R):
    """dipole along z with radial-spoke end plates ('top hats').
    zgap=None: one wire, returns its tag (feed at centre segment).
    zgap=d   : two arms ending at z=-d and z=+d (for the series model)."""
    tags = {}
    if zgap is None:
        tags["dip"] = st.wire([0, 0, -h / 2], [0, 0, h / 2], nd, "dipole", rad)
    else:
        tags["lo"] = st.wire([0, 0, -h / 2], [0, 0, -zgap], nd, "dipole", rad)
        tags["up"] = st.wire([0, 0, zgap], [0, 0, h / 2], nd, "dipole", rad)
    if r_p > 0:
        for s in (-1, 1):
            for i in range(n_sp):
                p = 2 * np.pi * (i + 0.5) / n_sp
                st.wire([0, 0, s * h / 2], [r_p * np.cos(p), r_p * np.sin(p), s * h / 2], sp_seg, "plates", rad)
    return tags


def model_loops(h=H_D, feed="outer", n_l=16, nd=5, r_p=0.03, n_sp=6, sp_seg=3, sense=+1,
                R=R_T, a=A_L, N=N_L, rad=WIRE_R, with_dipole=True, with_loops=True):
    """feed: 'outer' | 'inner' | 'top' | 'dual' (outer+inner) | 'quad'.
    Port groups (lists of (tag, seg)): 'dipole', 'loops' (+ 'loops2' for dual/quad).
    Positive loop current = +z on the outer side (sense=+1)."""
    st = Structure()
    ports = {}
    if with_dipole:
        t = add_dipole(st, h, nd, r_p, n_sp, sp_seg, rad=rad)
        ports["dipole"] = [(t["dip"], nd // 2 + 1)]
    if with_loops:
        ae = a * _eq_area(n_l)
        fp = dict(outer=[0], inner=[n_l // 2], top=[n_l // 4], dual=[0, n_l // 2],
                  quad=[0, n_l // 4, n_l // 2, 3 * n_l // 4])[feed]
        ports["loops"] = []
        if feed in ("dual", "quad"):
            ports["loops2"] = []
        for i in range(N):
            ph = 2 * np.pi * i / N
            psi = sense * (np.arange(n_l + 1) - 0.5) * 2 * np.pi / n_l
            rho = R + ae * np.cos(psi)
            pts = np.stack([rho * np.cos(ph), rho * np.sin(ph), ae * np.sin(psi)], -1)
            tg = st.polyline(pts, "coil", rad)
            for j in fp:
                key = "loops" if (j == 0 or feed not in ("dual", "quad")) else "loops2"
                ports[key].append((tg[j], 1))
    st.ports = ports
    return st


def model_series(h=H_D, n_t=16, nd=1, r_p=0.03, n_sp=6, sp_seg=3, sense=+1, adv=+1,
                 R=R_T, a=A_L, N=N_L, rad=WIRE_R, dz=0.006, R_ret=None, n_ret=48,
                 lead_seg=None, with_return=True,
                 caps_per_turn=0, C_turn=0.0, caps_ring=0, C_ring=0.0):
    """Series antenna, current path (for dipole current along +z):

    bottom plate -> lower arm -> lead L1 (z=-dz, axis -> rho=R-a)
      -> toroidal helix, N turns, starting on the inner equator (psi=pi),
         one full turn in phi (direction adv)
      -> jumper J1 (rho: R-a -> R_ret, z ~ 0)
      -> return ring, radius R_ret, one full turn in the direction -adv
      -> lead L2 (z=+dz, rho=R_ret -> axis) -> upper arm -> top plate.

    L1/L2 and J1/L2 are close antiparallel pairs (non-radiating).  The
    generator is on lead L1.  sense=+1: turn current is +z on the outer
    side when the dipole current is +z.
    """
    st = Structure()
    if R_ret is None:
        R_ret = np.sqrt(R ** 2 + a ** 2 / 2)
    ae = a * _eq_area(n_t)
    t = add_dipole(st, h, nd, r_p, n_sp, sp_seg, zgap=dz, rad=rad)
    seg_len = 2 * np.pi * ae / n_t
    ls = lead_seg or max(1, int(round((R - a) / seg_len)))
    # L1: axis -> inner equator at z=-dz
    tL1 = st.wire([0, 0, -dz], [R - ae, 0, -dz], ls, "leads", rad)
    # helix: psi from pi, N turns; z shifted linearly from -dz to 0
    n = N * n_t
    tt = np.arange(n + 1) / n
    psi = np.pi + sense * 2 * np.pi * N * tt
    ph = adv * 2 * np.pi * tt
    rho = R + ae * np.cos(psi)
    z = ae * np.sin(psi) - dz * (1 - tt)
    pts = np.stack([rho * np.cos(ph), rho * np.sin(ph), z], -1)
    pts[0] = [R - ae, 0, -dz]; pts[-1] = [R - ae, 0, 0]
    st.helix_tags = st.polyline(pts, "coil", rad)
    st.cap_turn, st.cap_ring = [], []
    if caps_per_turn:
        step = n_t // caps_per_turn
        st.cap_turn = [st.helix_tags[i] for i in range(step // 2, n, step)]
        for tg in st.cap_turn:
            st.load(tg, 1, C=C_turn)
    if with_return:
        nj = max(1, int(round((R_ret - (R - ae)) / seg_len)))
        st.wire([R - ae, 0, 0], [R_ret, 0, 0], nj, "leads", rad)                 # J1
        Re = R_ret * _eq_area(n_ret)
        # polygon with a vertex at phi=0 (start, z=0) and the same vertex at the end (z=+dz)
        tr = np.arange(n_ret + 1) / n_ret
        pr = -adv * 2 * np.pi * tr
        rr = np.full(n_ret + 1, Re); rr[0] = rr[-1] = R_ret
        ptsr = np.stack([rr * np.cos(pr), rr * np.sin(pr), dz * tr], -1)
        st.ret_tags = st.polyline(ptsr, "return", rad)
        if caps_ring:
            step = n_ret // caps_ring
            st.cap_ring = [st.ret_tags[i] for i in range(step // 2, n_ret, step)]
            for tg in st.cap_ring:
                st.load(tg, 1, C=C_ring)
        nl2 = max(1, int(round(R_ret / seg_len)))
        st.wire([R_ret, 0, dz], [0, 0, dz], nl2, "leads", rad)                    # L2
    else:
        # no return ring: come straight back from the helix end
        st.wire([R - ae, 0, 0], [R - ae, 0, dz], 1, "leads", rad)
        st.wire([R - ae, 0, dz], [0, 0, dz], ls, "leads", rad)
    st.feed = (tL1, 2)          # generator on lead L1, second segment from the axis
    st.source(*st.feed, 1.0)
    st.dip_ref = (t["lo"], nd // 2 + 1)
    st.ports = {"feed": [st.feed]}
    return st


# -------------------------------------------------------------------------
def solve_ports(st, f, targets, lossy=False):
    """Find the voltages of each port group so that the current of the FIRST
    port of each group equals targets[group]; return (Result, V dict, Y)."""
    groups = list(st.ports)
    Y = np.zeros((len(groups), len(groups)), complex)
    for j, g in enumerate(groups):
        r = st.solve(f, lossy, sources=[(t, s, 1.0 + 0j) for t, s in st.ports[g]])
        for i, gi in enumerate(groups):
            Y[i, j] = r.current(*st.ports[gi][0])
    V = np.linalg.solve(Y, np.array([targets[g] for g in groups], complex))
    src = [(t, s, V[j]) for j, g in enumerate(groups) for t, s in st.ports[g]]
    st.sources = src
    return st.solve(f, lossy, sources=src), dict(zip(groups, V)), Y


def model_twogen(h=H_D, n_t=16, nd=5, r_p=0.03, n_sp=6, sp_seg=3, sense=+1, adv=+1,
                 R=R_T, a=A_L, N=N_L, rad=WIRE_R, d=0.003, R_ret=None, n_ret=48,
                 caps_per_turn=0, X_turn=0.0, caps_ring=0, X_ring=0.0, f_caps=None,
                 with_dipole=True, with_coil=True, n_feeds=1):
    """Two-generator antenna (feeding option 2 of the report).

    * dipole with end plates, fed at its centre (port 'dipole');
    * CLOSED coil circuit: toroidal helix (N turns, one turn in phi) ->
      jumper -> return ring (one turn back) -> jumper -> helix start,
      fed by its own generator on the outer side of the first turn
      (port 'coil'; n_feeds > 1: identical generators on equally spaced turns).  Optional series capacitors: caps_per_turn in every
      turn (total reactance X_turn ohm per turn at f_caps) and caps_ring in
      the ring (total X_ring ohm), to keep the current uniform.
    """
    from .core import F0
    st = Structure()
    st.ports = {}
    w = 2 * np.pi * (f_caps or F0)
    if with_dipole:
        t = add_dipole(st, h, nd, r_p, n_sp, sp_seg, rad=rad)
        st.ports["dipole"] = [(t["dip"], nd // 2 + 1)]
    if not with_coil:
        return st
    if R_ret is None:
        R_ret = np.sqrt(R ** 2 + a ** 2 / 2)
    ae = a * _eq_area(n_t)
    n = N * n_t
    tt = np.arange(n + 1) / n
    psi = np.pi + sense * 2 * np.pi * N * tt
    ph = adv * 2 * np.pi * tt
    rho = R + ae * np.cos(psi)
    z = ae * np.sin(psi) + d * (2 * tt - 1)
    pts = np.stack([rho * np.cos(ph), rho * np.sin(ph), z], -1)
    pts[0] = [R - ae, 0, -d]; pts[-1] = [R - ae, 0, d]
    st.helix_tags = st.polyline(pts, "coil", rad)
    seg_len = 2 * np.pi * ae / n_t
    nj = max(1, int(round((R_ret - (R - ae)) / seg_len)))
    st.wire([R - ae, 0, d], [R_ret, 0, d], nj, "leads", rad)
    Re = R_ret * _eq_area(n_ret)
    tr = np.arange(n_ret + 1) / n_ret
    pr = -adv * 2 * np.pi * tr
    rr = np.full(n_ret + 1, Re); rr[0] = rr[-1] = R_ret
    ptsr = np.stack([rr * np.cos(pr), rr * np.sin(pr), d * (1 - 2 * tr)], -1)
    st.ret_tags = st.polyline(ptsr, "return", rad)
    st.wire([R_ret, 0, -d], [R - ae, 0, -d], nj, "leads", rad)
    # generator: outer side of the first turn (psi = 2 pi  -> index n_t/2)
    fidx = [n_t // 2 + j * (N // n_feeds) * n_t for j in range(n_feeds)]
    st.ports["coil"] = [(st.helix_tags[i], 1) for i in fidx]
    st.cap_turn, st.cap_ring = [], []
    if caps_per_turn:
        step = n_t // caps_per_turn
        idx = [i for i in range(step // 2, n, step)]
        st.cap_turn = [st.helix_tags[i] for i in idx if i not in fidx]
        C = caps_per_turn / (w * X_turn)
        for tg in st.cap_turn:
            st.load(tg, 1, C=C)
        st.C_turn = C
    if caps_ring:
        step = n_ret // caps_ring
        st.cap_ring = [st.ret_tags[i] for i in range(step // 2, n_ret, step)]
        C = caps_ring / (w * X_ring)
        for tg in st.cap_ring:
            st.load(tg, 1, C=C)
        st.C_ring = C
    return st
