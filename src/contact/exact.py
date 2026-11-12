"""Outward rational bounds for square roots."""
from fractions import Fraction as F
from math import isqrt

def sqrt_bounds(q, bits=200):
    """rational lo, hi with lo <= sqrt(q) <= hi for rational q >= 0"""
    q = F(q)
    s = 1 << bits
    num = q.numerator * s * s
    den = q.denominator
    # sqrt(num/den) = sqrt(num*den)/den
    t = isqrt(num * den)
    lo = F(t, den * s)
    hi = F(t + 1, den * s)
    assert lo * lo <= q <= hi * hi
    return lo, hi
