"""k-part partition certificates for the three-contact family in R^5.

For a dyadic parameter box Q the exact enclosure hd_enclose.enclose(Q, n) gives a rational outer body
    Omega_Q = {|x|^2 <= r_Q^+} ∩ ∩_i {|x - A_i|^2 <= b_i + |A_i|^2} ∩ F_n,
F_n = {|x_j| <= 1/2, j = 4..n-1} ∩ {x_n <= 1/2} (the frame in the orthogonal complement of the three contacts),
containing every member of the family represented by Q.  A certificate is a rational power diagram with k
cells; each cell polytope F_n ∩ cell_j ∩ (tangent half-spaces, each checked exactly to contain its ball of
Omega_Q) must be bounded with all vertex pairs at squared distance < 1.  The dimension is set with set_dim()."""
import os, sys, json, time
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import hd_enclose as E
from exact import sqrt_bounds
from shapes import Ball
from certify import float_vertices


def far_pairs(V, thr, top=400):
    """certify.far_pairs with the same output (pairs (i, j, d2) with d2 >= thr, the top ones by distance),
    the pair search is restricted to rows whose maximum reaches thr"""
    s = (V * V).sum(1); out = []; mx = 0.0
    for c in range(0, len(V), 2000):
        D2 = s[c:c+2000, None] + s[None, :] - 2 * V[c:c+2000] @ V.T
        rm = D2.max(1); mx = max(mx, rm.max())
        rows = np.nonzero(rm >= thr)[0]
        if len(rows):
            sub = D2[rows]; I, J = np.nonzero(sub >= thr)
            d = sub[I, J]; I = rows[I]
            if len(d) > top:
                sel = np.argpartition(-d, top)[:top]; I, J, d = I[sel], J[sel], d[sel]
            out.append(np.c_[I + c, J, d])
    if not out: return np.zeros((0, 3)), np.sqrt(mx)
    P = np.vstack(out); P = P[np.argsort(-P[:, 2])][:top]
    return P, np.sqrt(mx)

N = 5; DEN_N = 10**9; DEN_B = 10**12; HALF = Fr(1, 2)
# MIN_SEP > 0: a new tangent is skipped when an existing tangent of the same ball has a unit normal within MIN_SEP
# (avoids near-degenerate vertices). CDD = False: no cddlib fallback in check_cell (a failed attempt can be retried;
# a cell that needs cddlib is then reported as failed, which is always safe).
MIN_SEP = 0.0; CDD = True


def set_dim(n):
    global N
    N = int(n)


def ceil_frac(x, den=DEN_B):
    x = Fr(x); return Fr(-((-x.numerator * den) // x.denominator), den)


def outer(box):
    """exact balls [(centre, radius^2)] and planes of Omega_Q"""
    m = E.enclose(tuple(map(tuple, box)), N)
    balls = [([Fr(0)] * N, Fr(m['r_upper']))]
    for A, b in zip(m['contacts'], m['intercepts']):
        A = [Fr(x) for x in A]
        balls.append((A, Fr(b) + sum(x * x for x in A)))
    planes = E.frame_planes(N)
    return balls, planes, m


def float_body(balls, planes):
    fb = [(np.array([float(x) for x in c]), float(sqrt_bounds(R2, 60)[1])) for c, R2 in balls]
    fp = [(np.array([float(x) for x in a]), float(b)) for a, b in planes]
    return fb, fp


def tangent_ok(a, b, ball):
    c, R2 = ball
    s = b - sum(x * y for x, y in zip(a, c)); aa = sum(x * x for x in a)
    return s >= 0 and s * s >= R2 * aa


def make_tangent(u, ball):
    c, R2 = ball
    u = np.asarray(u, float); u = u / np.linalg.norm(u)
    a = [Fr(int(round(x * DEN_N)), DEN_N) for x in u]
    aa = sum(x * x for x in a)
    b = ceil_frac(sum(x * y for x, y in zip(a, c)) + sqrt_bounds(R2 * aa, 120)[1])
    assert tangent_ok(a, b, ball)
    return a, b


def cell_halfspaces(Aq, bq, i):
    return [([Aq[j][t] - Aq[i][t] for t in range(N)], bq[i] - bq[j]) for j in range(len(Aq)) if j != i]


def build_cell(balls, planes, Aq, bq, i, seeds=30, max_rounds=80, top=60, warm=None, verbose=False):
    fb, _ = float_body(balls, planes)
    shapes = [Ball(c, R) for c, R in fb]
    base = [(np.array([float(x) for x in a]), float(b)) for a, b in planes + cell_halfspaces(Aq, bq, i)]
    rng = np.random.default_rng(i)
    tans = {}; dirs = {}
    def add(kind, u):
        un = u / np.linalg.norm(u)
        key = (kind, tuple(np.round(un, 7)))
        if key in tans: return 0
        if MIN_SEP > 0:
            D = dirs.get(kind)
            if D is not None and len(D) and np.min(np.linalg.norm(np.array(D) - un, axis=1)) < MIN_SEP: return 0
            dirs.setdefault(kind, []).append(un)
        tans[key] = make_tangent(u, balls[kind]); return 1
    if warm:
        for kind, u in warm: add(kind, np.asarray(u, float))
        seeds = min(seeds, 6)
    for kind in range(len(balls)):
        for u in rng.standard_normal((seeds, N)): add(kind, u)
    dmax = 9.9
    for rnd in range(max_rounds):
        H = base + [(np.array([float(x) for x in a]), float(b)) for a, b in tans.values()]
        Vx = float_vertices(H)
        if Vx is None:
            for kind in range(len(balls)):
                for u in rng.standard_normal((15, N)): add(kind, u)
            continue
        if len(Vx) == 0: return tans, 0.0, True
        P, dmax = far_pairs(Vx, 1.0, top)
        if verbose: print(f'   cell {i} round {rnd}: {len(Vx)} verts, {len(tans)} tangents, dmax {dmax:.6f}', flush=True)
        if len(P) == 0: return tans, dmax, True
        idx = np.unique(P[:, :2].astype(int).ravel()); added = 0
        for v in Vx[idx]:
            for kind, sh in enumerate(shapes):
                if not sh.inside(v[None], 1e-10)[0]: added += add(kind, v - sh.c)
        if added == 0: return tans, dmax, False
    return tans, dmax, False




PRUNE_TOLS = (1e-4, 1e-3, 3e-3, 1e-2, 3e-2)


def prune(cell, tol):
    """a SUBSET of the tangent list: one tangent per (ball, direction rounded to tol), the latest one wins.
    Used only when near-identical tangents make the vertex enumeration degenerate; a subset of valid
    half-spaces gives a larger polytope that still contains Omega_Q ∩ cell, so this is always sound."""
    keep = {}
    for a, b, kind in cell:
        u = np.array([float(x) for x in a]); u /= np.linalg.norm(u)
        keep[(kind, tuple(np.round(u / tol).astype(np.int64)))] = (a, b, kind)
    return list(keep.values())


def check_cell(balls, planes, Aq, bq, i, cell):
    """exact: returns (ok, float diameter bound, #vertices).
    Every tangent is checked to contain its ball. Then the polytope of all half-spaces is enumerated exactly
    (adjacency-closure certificate). If that certificate cannot be completed (near-degenerate vertices), the
    same is tried for pruned subsets of the tangents and finally cddlib is run on the full system."""
    from fastcheck import fast_vertices
    from check_certificate import exact_vertices, max_sq_dist
    base = list(planes) + cell_halfspaces(Aq, bq, i)
    for a, b, kind in cell:
        if not tangent_ok(a, b, balls[kind]): return False, None, 0
    H = base + [(a, b) for a, b, _ in cell]

    def judge(V):
        if not V: return True, 0.0, 0
        ex, fl = max_sq_dist(V)
        if ex is not None and ex >= 1: return False, float(ex) ** .5, len(V)
        return True, float(fl) ** .5, len(V)
    V = fast_vertices(H, N)
    if V is not None: return judge(V)
    for tol in PRUNE_TOLS:
        sub = prune(cell, tol)
        V = fast_vertices(base + [(a, b) for a, b, _ in sub], N)
        if V is not None:
            res = judge(V)
            if res[0]: return res
    if not CDD: return False, None, 0
    V = exact_vertices(H, N)
    if V is None: return False, None, 0
    return judge(V)
