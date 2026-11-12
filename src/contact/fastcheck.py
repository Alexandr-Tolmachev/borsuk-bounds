"""Exact vertex enumeration certified by adjacency closure (fast alternative to a cddlib GMP run).

Input: K = { x in R^n : G x <= g } with rational data. A floating-point enumeration (Qhull) proposes vertices
as index sets I of n tight constraints. Each proposal is then verified in exact arithmetic:

  (a) G_I is nonsingular (its determinant is nonzero modulo the prime 2^61 - 1) and x_I = G_I^{-1} g_I is computed
      exactly (fraction-free Bareiss elimination on integers), and G_I x_I = g_I is re-checked exactly;
  (b) G_j x_I < g_j strictly for every j not in I  (x_I is a vertex of K and is nondegenerate);
  (c) closure: for every verified I and every i in I there is another verified J != I with I \\ {i} ⊂ J.

Claim: if (a)-(c) hold for a nonempty family S, then {x_I : I in S} is the complete vertex set of K and K is
bounded (hence K = conv of these points and diam K = max distance between them).
Proof. Each x_I is a nondegenerate vertex, so its incident edges are exactly the n one-dimensional faces
F_i = K ∩ {G_{I\\{i}} x = g_{I\\{i}}}. By (c) F_i contains another vertex x_J of K; a vertex of K lying on the
line of F_i is an endpoint of the segment/ray F_i, so F_i = [x_I, x_J] is bounded and its other endpoint is in S.
K has a vertex, so it is pointed; the graph formed by the vertices and bounded edges of a pointed polyhedron
is connected, so a vertex set closed under bounded edges, starting from a vertex, contains all vertices. If K
were unbounded it would have an unbounded edge, and every edge of a pointed polyhedron contains a vertex;
the edges at vertices of S are all bounded, so K is bounded. ∎

Float arithmetic is used only to (i) propose the index sets and (ii) skip the exact test in (b) when the
floating-point slack exceeds a rigorous rounding-error bound; every other decision is exact. If anything fails
(degenerate or missed vertices, an empty or lower-dimensional K, Qhull errors) the caller falls back to cddlib.
"""
import numpy as np
from fractions import Fraction as Fr
from math import lcm
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection

U = 2.0 ** -53

def int_rows(H):
    """Scale every inequality a.x <= b to integer coefficients (row-wise lcm of denominators)."""
    out = []
    for a, b in H:
        d = 1
        for x in list(a) + [b]: d = lcm(d, Fr(x).denominator)
        out.append(([int(Fr(x) * d) for x in a], int(Fr(b) * d)))
    return out

def solve_int(M, rhs):
    """Exact solution of M x = rhs (integers) as (X, D) with x = X / D, D > 0; None if singular. Bareiss."""
    n = len(M)
    A = [list(M[r]) + [rhs[r]] for r in range(n)]
    prev = 1
    for k in range(n):
        p = next((r for r in range(k, n) if A[r][k] != 0), None)
        if p is None: return None
        if p != k: A[k], A[p] = A[p], A[k]
        akk = A[k][k]
        for r in range(k + 1, n):
            ark = A[r][k]
            Ar = A[r]; Ak = A[k]
            for c in range(k + 1, n + 1):
                Ar[c] = (akk * Ar[c] - ark * Ak[c]) // prev
            Ar[k] = 0
        prev = akk
    D = A[n - 1][n - 1]          # = det (up to the sign of the row swaps), all divisions above are exact
    # back substitution, keeping x = X / D with integer X
    X = [0] * n
    for k in range(n - 1, -1, -1):
        s = A[k][n] * D - sum(A[k][c] * X[c] for c in range(k + 1, n))
        q, r = divmod(s, A[k][k])
        if r: return None        # cannot happen for a nonsingular Bareiss form; be safe
        X[k] = q
    if D < 0: X = [-x for x in X]; D = -D
    return X, D

P61 = (1 << 61) - 1                                    # a Mersenne prime

def det_mod_p(M, p=P61):
    """Determinant of an integer matrix modulo p (Gaussian elimination in GF(p)). Nonzero => det != 0."""
    A = [[x % p for x in r] for r in M]; n = len(A); det = 1
    for k in range(n):
        r = next((r for r in range(k, n) if A[r][k]), None)
        if r is None: return 0
        if r != k: A[k], A[r] = A[r], A[k]; det = -det
        det = det * A[k][k] % p; inv = pow(A[k][k], p - 2, p)
        for r in range(k + 1, n):
            f = A[r][k] * inv % p
            if f:
                for c in range(k, n): A[r][c] = (A[r][c] - f * A[k][c]) % p
    return det % p

def fast_vertices(H, n, return_reason=False):
    """Returns the exact vertex list (Fractions) of {x : a.x <= b for (a, b) in H}, certified as above,
    or None when the certificate could not be completed (then use cddlib)."""
    def fail(r): return (None, r) if return_reason else None
    Hi = int_rows(H)
    G = [a for a, b in Hi]; g = [b for a, b in Hi]
    Gf = np.array([[float(Fr(x)) for x in a] for a, b in H]); gf = np.array([float(Fr(b)) for a, b in H])
    if not np.isfinite(Gf).all() or not np.isfinite(gf).all():
        return fail('nonfinite constraint data')
    nrm = np.linalg.norm(Gf, axis=1)
    res = linprog(np.r_[np.zeros(n), -1], A_ub=np.c_[Gf, nrm], b_ub=gf, bounds=[(-5, 5)] * n + [(0, 1)], method='highs')
    if res.status != 0 or res.x[-1] < 1e-7: return fail('no interior point')
    c = res.x[:n]
    try:
        hs = HalfspaceIntersection(np.c_[Gf, -gf], c)
    except Exception as e:
        return fail(f'qhull {e}')
    S = {}
    for f in hs.dual_facets:
        if len(f) != n: return fail('non-simplicial dual facet')
        S[frozenset(int(x) for x in f)] = None
    if not S:
        return fail('empty vertex proposal')
    m = len(H)
    absG = np.abs(Gf); absg = np.abs(gf)
    for I in S:
        Il = sorted(I)
        sol = solve_int([G[i] for i in Il], [g[i] for i in Il])
        if sol is None: return fail('singular')
        X, D = sol
        # independent re-checks of the solver: G_I x = g_I exactly, and det G_I != 0 (nonzero modulo a prime)
        if any(sum(G[i][t] * X[t] for t in range(n)) != g[i] * D for i in Il): return fail('solver residual')
        if det_mod_p([G[i] for i in Il]) == 0: return fail('det = 0 mod p')
        xf = np.array([x / D for x in X])                 # float(X/D) via true division of ints: correctly rounded
        slack = gf - Gf @ xf
        tol = 1e-12 * (absg + absG @ np.abs(xf) + 1.0)    # >> (n+3)·2u·(|g|+|G||x|): rigorous screen
        if not np.isfinite(xf).all() or not np.isfinite(slack).all() or not np.isfinite(tol).all():
            return fail('nonfinite inactivity filter')
        mask = np.ones(m, bool); mask[Il] = False
        for j in np.nonzero(mask & (slack <= tol))[0]:    # exact test only where the float screen is inconclusive
            if g[j] * D - sum(G[j][t] * X[t] for t in range(n)) <= 0: return fail('infeasible or degenerate vertex')
        S[I] = (X, D)
    # closure
    cnt = {}
    for I in S:
        for i in I:
            k = I - {i}; cnt[k] = cnt.get(k, 0) + 1
    if any(v != 2 for v in cnt.values()): return fail('closure')
    V = [[Fr(x, D) for x in X] for X, D in S.values()]
    return (V, 'ok') if return_reason else V
