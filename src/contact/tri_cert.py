"""Six-part partition certificates for the three-contact family in R^4.

For a dyadic parameter box Q (coordinates as in Arman-Bondarenko-Prymak-Radchenko, r in [1/4, 2/5]) the
exact enclosure six_balls.enclose(Q) gives a rational outer body
    Omega_Q = {|x|^2 <= r_Q^+} ∩ ∩_i {|x - A_i|^2 <= b_i + |A_i|^2} ∩ {x4 <= 1/2}
containing every member of the family represented by Q (the half-space x4 <= 1/2 comes from reflecting the
set in the hyperplane of the three contacts and the circumcentre).  A certificate is a rational power diagram
(a_j, beta_j), j = 1..6, and for every cell a list of rational half-spaces, each tagged with the ball of
Omega_Q it contains.  The exact check: every tagged half-space contains its ball; the polytope
{x4 <= 1/2} ∩ cell_j ∩ (half-spaces) is bounded and all its vertex pairs are at squared distance < 1.
Then diam(Omega_Q ∩ cell_j) < 1 for every j, so every set of diameter 1 represented by Q splits into six
parts of smaller diameter."""
import os, sys, json, time
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import six_balls as S
from exact import sqrt_bounds
from shapes import Ball
from certify import float_vertices, far_pairs

N = 4; DEN_N = 10**9; DEN_B = 10**12; HALF = Fr(1, 2)


def ceil_frac(x, den=DEN_B):
    x = Fr(x); return Fr(-((-x.numerator * den) // x.denominator), den)


def outer(box):
    """exact balls [(centre, radius^2)] and planes of Omega_Q"""
    m = S.enclose(tuple(map(tuple, box)))
    balls = [([Fr(0)] * N, Fr(m['r_upper']))]
    for A, b in zip(m['contacts'], m['intercepts']):
        A = [Fr(x) for x in A]
        balls.append((A, Fr(b) + sum(x * x for x in A)))
    planes = [([Fr(0), Fr(0), Fr(0), Fr(1)], HALF)]
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
    tans = {}
    def add(kind, u):
        key = (kind, tuple(np.round(u / np.linalg.norm(u), 7)))
        if key in tans: return 0
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




PRUNE_TOLS = (1e-4, 1e-3, 3e-3)


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
    V = exact_vertices(H, N)
    if V is None: return False, None, 0
    return judge(V)
