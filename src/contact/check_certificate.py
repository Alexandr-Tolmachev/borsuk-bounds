"""Rational vertex enumeration and certified squared-distance comparisons."""
from fractions import Fraction as F
import numpy as np
import cdd.gmp as g

def exact_vertices(H, n):
    rows = [[b] + [-x for x in a] for a, b in H]
    mat = g.matrix_from_array(rows, rep_type=g.RepType.INEQUALITY)
    poly = g.polyhedron_from_matrix(mat)
    gen = g.copy_generators(poly)
    V = []
    for row in gen.array:
        if row[0] != 1:
            return None          # a ray: unbounded (or lineality)
        V.append([F(x) for x in row[1:]])
    if gen.lin_set:
        return None
    return V

def max_sq_dist(V):
    """Exact test of max squared distance among vertices V (rationals).
    Float screen: coordinates are rounded to doubles (|x_t| < 4, relative error <= u = 2^-53) and the squared
    distance of a pair is computed as |x|^2 + |y|^2 - 2 x.y. Each of the three dot products has absolute error
    <= gamma_n * 16n (gamma_n = nu/(1-nu)), the input rounding adds <= 4u * 16n, the two final operations
    <= 2u * 64n; for n <= 8 the total is < 1e-11. Pairs whose float value exceeds 1 - 1e-9 are re-checked
    exactly, so a pair with exact squared distance >= 1 cannot be missed. Returns (exact max over the
    re-checked pairs or None, float max)."""
    Vf = np.array([[float(x) for x in v] for v in V])
    assert Vf.shape[1] <= 8 and np.abs(Vf).max() < 4
    sq = (Vf * Vf).sum(1)
    fmax = 0.0; cand = []
    B = 1024
    for c in range(0, len(Vf), B):              # block rows c..c+B against columns c.. (pairs i < j)
        D2 = Vf[c:c + B] @ Vf[c:].T
        D2 *= -2.0; D2 += sq[None, c:]; D2 += sq[c:c + B, None]
        fmax = max(fmax, D2.max())
        I, J = np.nonzero(D2 > 1 - 1e-9)
        cand += [(i + c, j + c) for i, j in zip(I, J) if i < j]
    best = None
    for i, j in cand:
        d = sum((a - b) ** 2 for a, b in zip(V[i], V[j]))
        if best is None or d > best: best = d
    return best, fmax
