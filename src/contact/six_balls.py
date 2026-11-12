"""Rational interval enclosures for the four-dimensional contact family.

The interval arithmetic is derived from seven_balls_single.py by
Andrii Arman, Andriy Bondarenko, Andriy Prymak and Danylo Radchenko:
https://github.com/andriyprm/4d7balls
"""
from fractions import Fraction as F
from math import isqrt

R_LO, R_HI = F(1, 4), F(2, 5)
DELTA_LO, DELTA_HI = F(1, 3)-R_LO, R_HI-F(1, 3)

def require(condition, message):
    """Never use Python assert for certificate acceptance (python -O removes it)."""
    if not condition:
        raise ValueError(message)

def ceil_div(a, b):
    """Integer ceiling of a/b; used below only with b > 0."""
    return -((-a) // b)

class Interval:
    """Closed interval with integer endpoints on the grid 2**(-44) Z.

    All intermediate operations are rounded OUTWARD. For multiplication and
    division, the extrema occur among the four endpoint pairs. Division is
    allowed only when the denominator interval excludes zero. The method
    positive() encloses max(0,x); it does not merely discard negative values.
    """
    SCALE = 1 << 44
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        a = F(lo)
        b = a if hi is None else F(hi)
        self.lo = a.numerator * self.SCALE // a.denominator
        self.hi = ceil_div(b.numerator * self.SCALE, b.denominator)
        require(self.lo <= self.hi, 'Reversed interval')

    @classmethod
    def raw(cls, lo, hi):
        require(type(lo) is int and type(hi) is int and lo <= hi,
                'Invalid fixed-point interval')
        obj = object.__new__(cls)
        obj.lo, obj.hi = lo, hi
        return obj

    @classmethod
    def coerce(cls, x):
        return x if isinstance(x, cls) else cls(x)

    def __add__(self, other):
        other = self.coerce(other)
        return self.raw(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__

    def __neg__(self):
        return self.raw(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) + (-self)

    def __mul__(self, other):
        other = self.coerce(other)
        values = [a*b for a in (self.lo, self.hi)
                  for b in (other.lo, other.hi)]
        return self.raw(min(values)//self.SCALE,
                        ceil_div(max(values), self.SCALE))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError('Interval denominator contains zero')
        if other.hi < 0:
            return (-self)/(-other)
        lows = [a*self.SCALE//b for a in (self.lo, self.hi)
                for b in (other.lo, other.hi)]
        highs = [ceil_div(a*self.SCALE, b) for a in (self.lo, self.hi)
                 for b in (other.lo, other.hi)]
        return self.raw(min(lows), max(highs))

    def __rtruediv__(self, other):
        return self.coerce(other)/self

    def square(self):
        lo = 0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        return self.raw(lo//self.SCALE,
                        ceil_div(max(self.lo**2, self.hi**2), self.SCALE))

    def sqrt(self):
        require(self.lo >= 0, 'Negative square-root argument')
        lo = isqrt(self.lo*self.SCALE)
        square = self.hi*self.SCALE
        hi = isqrt(square)
        return self.raw(lo, hi if hi*hi == square else hi+1)

    def positive(self):
        return self.raw(max(0, self.lo), max(0, self.hi))

    def intersect(self, other):
        other = self.coerce(other)
        return self.raw(max(self.lo, other.lo), min(self.hi, other.hi))

    def left(self):
        return self.raw(self.lo, self.lo)

    def right(self):
        return self.raw(self.hi, self.hi)

    def endpoints(self):
        return F(self.lo, self.SCALE), F(self.hi, self.SCALE)

def axis_bounds(axis):
    k, level = axis
    require(type(k) is int and type(level) is int
            and 0 <= level <= 128 and 0 <= k < 2**level,
            'Invalid dyadic coordinate')
    return F(k, 2**level), F(k+1, 2**level)

def parameter_r_interval(axis):
    """r=1/3+t|t|, t=x sqrt(DELTA_HI)-(1-x)sqrt(DELTA_LO).

    t and r increase with x. Evaluating the two endpoints with outward
    intervals therefore bounds the entire interval, including t=0.
    The identity t|t|=t_+**2-(-t)_+**2 needs no division or sign decision.
    """
    I = Interval
    low, high = map(I, axis_bounds(axis))
    dl, du = I(DELTA_LO).sqrt(), I(DELTA_HI).sqrt()

    def at(x):
        t = x*du-(1-x)*dl
        r = (I(F(1, 3))+t.positive().square()-(-t).positive().square())
        return t, r.intersect(I(R_LO, R_HI))

    tl, rl = at(low)
    tu, ru = at(high)
    return I.raw(rl.lo, ru.hi), (tl, rl), (tu, ru)

def rational_outer_model(r_interval, contact_intervals):
    """Enclose uncertain contacts by rational centres and outward L2 errors.

    With centre A_i, coordinate half-widths e_ij and r <= r_upper,
    |(P_i-A_i).X| <= sqrt(r_upper * sum_j e_ij**2).
    """
    contacts = [[(lo+hi)/2 for lo, hi in row] for row in contact_intervals]
    rl, ru = r_interval
    # Every point X satisfies |X|^2 <= r <= ru, and |P_i - A_i| <= |e_i|_2,
    # so |(P_i - A_i).X| <= sqrt(ru * sum_j e_ij^2) =: E_i (rational upper bound, outward rounding).
    def sqrt_upper(q, bits=40):
        q = F(q); scale = 1 << bits
        m = isqrt(q.numerator*q.denominator*scale*scale)
        return F(m+1, q.denominator*scale)
    errors = [sqrt_upper(ru*sum(((hi-lo)/2)**2 for lo, hi in row)) for row in contact_intervals]
    require(all(e*e >= ru*sum(((hi-lo)/2)**2 for lo, hi in row)
                for e, row in zip(errors, contact_intervals)), 'sqrt upper bound failed')
    return {'r_upper': ru, 'contacts': contacts,
            'intercepts': [1-rl+2*e for e in errors],
            'uncertainty': ru-rl+max(errors),
            'halfwidths': [[(hi-lo)/2 for lo, hi in row] for row in contact_intervals], 'r_lower': rl}

def enclose(box):
    """One exact enclosure for the WHOLE computational interval.

    With (x,y,z,w) in [0,1]^4, set
      t=x sqrt(r_hi-1/3)-(1-x)sqrt(1/3-r_lo), r=1/3+t|t|,
      h=(1-y)sqrt(r/6)+y max(t,0), tau=r-h*h,
      a0=sqrt(max(0,1-1/(4tau))), a=(1-z)/2+z*a0,
      b=(1-w)sqrt((1-a)/2)+w*a.
    Then 5r/6 <= tau <= min(r,1/3). The contact coordinates are
      (sqrt(tau),0,h,0),
      (sqrt(tau)(2a*a-1), 2sqrt(tau)a sqrt(1-a*a),h,0),
      (sqrt(tau)(2b*b-1),-2sqrt(tau)b sqrt(1-b*b),h,0).

    Height increases with x and decreases with y: both height endpoints
    increase with x, and max(t,0) <= sqrt(r/6) since r <= 2/5. Thus corner
    bounds avoid unnecessary dependency loss in the height interpolation.
    All remaining intersections below use proved analytic bounds, not
    numerical observations. Fixed-point roundoff is always OUTWARD.
    """
    require(len(box) == 4, 'A parameter box must have four coordinates')
    I = Interval
    r, (tl, rl), (tu, ru) = parameter_r_interval(box[0])
    y, z, w = [I(*axis_bounds(axis)) for axis in box[1:]]
    hl = (1-y.right())*(rl/6).sqrt()+y.right()*tl.positive()
    hu = (1-y.left())*(ru/6).sqrt()+y.left()*tu.positive()
    h = I.raw(hl.lo, hu.hi).positive()
    tau = (r-h.square()).intersect(I(F(5, 24), F(1, 3)))
    tau = I.raw(max(tau.lo, (r*F(5, 6)).lo), min(tau.hi, r.hi))
    a0 = (1-1/(4*tau)).positive().sqrt().intersect(I(0, F(1, 2)))
    # For each actual tau, a decreases with z. a0 increases with tau.
    a = I.raw(((1-z.right())*F(1, 2)+z.right()*a0.left()).lo,
              ((1-z.left())*F(1, 2)+z.left()*a0.right()).hi)
    a = a.intersect(I(0, F(1, 2)))
    base = ((1-a)/2).sqrt()
    # base >= a for 0 <= a <= 1/2, so b decreases with w. Do NOT assume
    # b is monotone in a: independent endpoint bounds are used here.
    b = I.raw(((1-w.right())*base.left()+w.right()*a.left()).lo,
              ((1-w.left())*base.right()+w.left()*a.right()).hi)
    b = b.intersect(I(0, 1))

    def direction(c):
        # Both a and b lie in [0,1/sqrt(2)]. On that interval, 2c*c-1
        # and 2c sqrt(1-c*c) increase. Evaluate endpoint squares to avoid
        # treating c and sqrt(1-c*c) as unrelated interval variables.
        z = c.square().intersect(I(0, F(1, 2)))
        sl = (4*z.left()*(1-z.left())).positive().sqrt()
        su = (4*z.right()*(1-z.right())).positive().sqrt()
        return 2*z-1, I.raw(sl.lo, su.hi)

    ca, sa = direction(a)
    cb, sb = direction(b)
    v, zero = tau.sqrt(), I(0)
    raw = [[v, zero, h, zero], [v*ca, v*sa, h, zero],
           [v*cb, -v*sb, h, zero]]
    return rational_outer_model(r.endpoints(),
                               [[value.endpoints() for value in row] for row in raw])
