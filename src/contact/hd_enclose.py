"""Rational enclosures and frame inequalities for the contact family.

Uses the interval arithmetic of Arman, Bondarenko, Prymak and Radchenko
through six_balls. The contact radius satisfies c*r <= tau <= min(r,1/3),
where c=2*(n+1)/(3*n); see the three-contact normal form in the article.
"""
from fractions import Fraction as F
import six_balls as S
I = S.Interval

def constants(n):
    S.require(n in (4, 5, 6), 'dimension must be 4, 5 or 6')
    c = F(2 * (n + 1), 3 * n)
    r_lo, r_hi = F(1, 4), F(n, 2 * (n + 1))
    S.require(c * r_hi == F(1, 3), 'r_hi = 1/(3c)')
    return c, r_lo, r_hi

def parameter_r_interval(axis, n):
    c, r_lo, r_hi = constants(n)
    low, high = map(I, S.axis_bounds(axis))
    dl, du = I(F(1, 3) - r_lo).sqrt(), I(r_hi - F(1, 3)).sqrt()

    def at(x):
        t = x * du - (1 - x) * dl
        r = I(F(1, 3)) + t.positive().square() - (-t).positive().square()
        return t, r.intersect(I(r_lo, r_hi))
    tl, rl = at(low)
    tu, ru = at(high)
    return I.raw(rl.lo, ru.hi), (tl, rl), (tu, ru)

def contact_intervals(box, n):
    """interval enclosures (as Interval objects) of r and of the n coordinates of P1, P2, P3 over the box"""
    S.require(len(box) == 4, 'A parameter box must have four coordinates')
    c, r_lo, r_hi = constants(n)
    # (1 - c) r = (n-2) r / (3n); written with integer operations so that n = 4 reproduces r/6 exactly
    r, (tl, rl), (tu, ru) = parameter_r_interval(box[0], n)
    y, z, w = [I(*S.axis_bounds(axis)) for axis in box[1:]]
    # h decreases with y and both of its endpoints increase with x (max(t,0) <= sqrt((1-c) r), see docstring)
    hl = (1 - y.right()) * ((rl * (n - 2)) / (3 * n)).sqrt() + y.right() * tl.positive()
    hu = (1 - y.left()) * ((ru * (n - 2)) / (3 * n)).sqrt() + y.left() * tu.positive()
    h = I.raw(hl.lo, hu.hi).positive()
    tau = (r - h.square()).intersect(I(c * r_lo, F(1, 3)))
    tau = I.raw(max(tau.lo, (r * c).lo), min(tau.hi, r.hi))
    a0 = (1 - 1 / (4 * tau)).positive().sqrt().intersect(I(0, F(1, 2)))
    a = I.raw(((1 - z.right()) * F(1, 2) + z.right() * a0.left()).lo,
              ((1 - z.left()) * F(1, 2) + z.left() * a0.right()).hi)
    a = a.intersect(I(0, F(1, 2)))
    base = ((1 - a) / 2).sqrt()
    b = I.raw(((1 - w.right()) * base.left() + w.right() * a.left()).lo,
              ((1 - w.left()) * base.right() + w.left() * a.right()).hi)
    b = b.intersect(I(0, 1))

    def direction(cc):
        zz = cc.square().intersect(I(0, F(1, 2)))
        sl = (4 * zz.left() * (1 - zz.left())).positive().sqrt()
        su = (4 * zz.right() * (1 - zz.right())).positive().sqrt()
        return 2 * zz - 1, I.raw(sl.lo, su.hi)
    ca, sa = direction(a)
    cb, sb = direction(b)
    v, zero = tau.sqrt(), I(0)
    pad = [zero] * (n - 3)
    raw = [[v, zero, h] + pad, [v * ca, v * sa, h] + pad, [v * cb, -v * sb, h] + pad]
    return r, raw

def enclose(box, n):
    """rational outer model (L2 error bound): r_upper, contacts A_i, intercepts 1 - r_lower + 2 E_i, ..."""
    r, raw = contact_intervals(box, n)
    return S.rational_outer_model(r.endpoints(), [[value.endpoints() for value in row] for row in raw])

def frame_planes(n):
    """|x_j| <= 1/2 for 0-based j = 3..n-2 and x_{n-1} <= 1/2 (rational)"""
    H = []
    for j in range(3, n - 1):
        e = [F(0)] * n; e[j] = F(1); H.append((e, F(1, 2)))
        e = [F(0)] * n; e[j] = F(-1); H.append((e, F(1, 2)))
    e = [F(0)] * n; e[n - 1] = F(1); H.append((e, F(1, 2)))
    return H
