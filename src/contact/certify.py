"""Numerical proposals for tangent refinement; acceptance is checked exactly."""
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection

def float_vertices(H):
    """H: list of (a,b) a.x<=b (floats). returns vertex array or None if unbounded/empty"""
    Aa = np.array([a for a, b in H], float); bb = np.array([b for a, b in H], float)
    n = Aa.shape[1]
    nrm = np.linalg.norm(Aa, axis=1)
    res = linprog(np.r_[np.zeros(n), -1], A_ub=np.c_[Aa, nrm], b_ub=bb,
                  bounds=[(-5, 5)] * n + [(0, 1)], method='highs')
    if res.status == 2: return np.zeros((0, n))   # infeasible: the (float) polytope is empty
    if res.status != 0: return None
    if res.x[-1] < 1e-9: return np.zeros((0, n))
    c = res.x[:n]
    if np.any(np.abs(c) > 4.9): return None
    hs = HalfspaceIntersection(np.c_[Aa, -bb], c)
    V = hs.intersections
    if np.abs(V).max() > 4: return None
    return np.unique(np.round(V, 12), axis=0)

def far_pairs(V, thr, top=400):
    """pairs (i,j,d2) with squared distance >= thr, the top ones by distance, chunked"""
    s = (V * V).sum(1); out = []; mx = 0.0
    for c in range(0, len(V), 2000):
        D2 = s[c:c+2000, None] + s[None, :] - 2 * V[c:c+2000] @ V.T
        mx = max(mx, D2.max())
        I, J = np.nonzero(D2 >= thr)
        if len(I):
            d = D2[I, J]
            if len(d) > top:
                sel = np.argpartition(-d, top)[:top]; I, J, d = I[sel], J[sel], d[sel]
            out.append(np.c_[I + c, J, d])
    if not out: return np.zeros((0, 3)), np.sqrt(mx)
    P = np.vstack(out); P = P[np.argsort(-P[:, 2])][:top]
    return P, np.sqrt(mx)
