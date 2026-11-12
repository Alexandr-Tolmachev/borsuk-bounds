"""Exact arithmetic supporting the analytic proofs; standard library only.

Checks identities, explicit witnesses, admissibility and numerical bounds.
The all-dimensional geometric proofs are in the article, not inferred from tests.
No sampling, optimization, Qhull or contact-certificate replay is used.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import re
import argparse

BLOCKS = [(1,2,F(0)), (5,30,F(3,40)), (6,44,F(1,19)), (7,88,F(0)),
          (8,156,F(0)), (9,288,F(1,24)), (11,1244,F(1,35)), (12,2312,F(4225,125329))]
TARGETS = {4:6,5:12,6:31,7:45,8:89,9:157,10:289,11:625,12:1245,13:2313,
           14:4977,15:7745,16:13729,17:24337,18:48673,19:97345}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def slack(n, mu):
    p, q = mu.numerator, mu.denominator
    return 2*n*q*q-(n-1)*(q+(n+1)*p)**2


def best(n):
    usable = [(m, count) for m, count, mu in BLOCKS if slack(n, mu)>0]
    dp, witnesses = [1], [()]
    for dimension in range(1,n):
        value, witness = min((dp[dimension-m]*count, witnesses[dimension-m]+(m,))
                             for m,count in usable if m<=dimension)
        dp.append(value)
        witnesses.append(witness)
    return dp[-1]+1, witnesses[-1]


def root_parameters(m, weight):
    ell = (m-3-weight)/(weight-1)
    return ((weight-2)**2/(1+(weight-2)**2), max(F(0),(2*ell-m+2)/(ell**2+m-1)))


def star_configuration(m):
    vectors = [e for e in product((-1,1), repeat=m) if e.count(-1)%2==0]
    for i in range(m-1):
        for a,b in product((-2,2), repeat=2):
            v = [0]*m
            v[i],v[-1] = a,b
            vectors.append(tuple(v))
    return vectors


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def in_region(x,v,vectors):
    return all(dot(x,v)>=dot(x,w) for w in vectors)


class Polynomial:
    """Sparse integer polynomials in six variables, for two proof identities."""
    def __init__(self, terms):
        self.terms = {key:value for key,value in terms.items() if value}

    @staticmethod
    def lift(value):
        return value if isinstance(value,Polynomial) else Polynomial({(0,)*6:value})

    def __add__(self, other):
        result = self.terms.copy()
        for key,value in self.lift(other).terms.items():
            result[key] = result.get(key,0)+value
        return Polynomial(result)

    __radd__ = __add__

    def __mul__(self, other):
        result = {}
        for a,x in self.terms.items():
            for b,y in self.lift(other).terms.items():
                key = tuple(i+j for i,j in zip(a,b))
                result[key] = result.get(key,0)+x*y
        return Polynomial(result)

    __rmul__ = __mul__

    def __sub__(self, other):
        return self+(-1)*self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+(-1)*self

    def __pow__(self, exponent):
        result = self.lift(1)
        for _ in range(exponent):
            result = result*self
        return result


def check_star_algebra():
    k,d,s,t,u,v = [Polynomial({tuple(int(i==j) for i in range(6)):1}) for j in range(6)]
    # s,t = sigma,sigma'; u,v = delta,delta'.
    lhs = (s+(k-3)*u)*(t+(k-3)*v)-((d-1)*u*v+v*s+u*t+s*t)
    rhs = (k-4)*(u*t+v*s)+((k-5)*(k-2)+(k-d))*u*v
    require(not (lhs-rhs).terms, 'star nonnegative-sigma identity')
    # Here s,t stand for the original last coordinates.
    lhs = s*t-u*t-v*s+(2*k-5-d)*u*v
    rhs = (s-u)*(t-v)+(2*k-6-d)*u*v
    require(not (lhs-rhs).terms, 'star negative-sigma identity')
    require(not ((2*k-6-(k-1))-(k-5)).terms, 'star parity case d=k-1')
    require(not ((2*k-6-k)-(k-6)).terms, 'star parity case d=k')
    print('PASS: symbolic star identities; nonnegativity uses k>=5 and parity as in the proof')


def check_table(article=None):
    if article:
        lines = (article/'sections/introduction.tex').read_text().splitlines()
        printed = {int(t.split('&')[0]):int(t.split('&')[2]) for t in lines if re.match(r'^\d+&',t)}
        require(printed==TARGETS, 'Table 1 differs from the checked targets')
        print('PASS: Table 1 agrees with the certificate targets')
    for n in range(6,20):
        value,witness = best(n)
        require(value==TARGETS[n], f'block bound in dimension {n}')
        print(f'PASS: b({n}) <= {value}; block dimensions {witness}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--article', type=Path, help='also compare with Table 1 in the separate LaTeX sources')
    args = parser.parse_args()
    print('Exact arithmetic supporting the analytic proofs.')
    check_star_algebra()
    for m,count,mu in [(7,88,F(0)), (8,156,F(0)), (9,288,F(1,24))]:
        require(max(F(0),F(m-8,(m-5)**2+m-1))==mu, 'star tolerance')
        vectors = star_configuration(m)
        require(len(vectors)==len(set(vectors))==count, 'star cardinality')
        root = (2,)+(0,)*(m-2)+(2,)
        x,y = (1,1)+(0,)*(m-2),(1,-1)+(0,)*(m-2)
        require(in_region(x,root,vectors) and in_region(y,root,vectors) and dot(x,y)==0, 'zero witness')
        ell = m-5
        x,y = (-1,ell)+(1,)*(m-2),(ell,-1)+(1,)*(m-2)
        require(in_region(x,(1,)*m,vectors) and in_region(y,(1,)*m,vectors), 'sign-region witnesses')
        require(F(dot(x,y),dot(x,x))==F(m-2-2*ell,ell*ell+m-1), 'witness value')
        print(f'PASS: star block ({m},{count},{mu}), explicit boundary witnesses')
    x,y,v = (1,1,1,1,-1,1),(1,-1,-1,-1,1,1),(2,0,0,0,0,2)
    require(in_region(x,v,star_configuration(6)) and in_region(y,v,star_configuration(6))
            and F(dot(x,y),dot(x,x))==F(-1,3), 'm=6 excluded example')
    require(F(1,3)**2/(2+F(1,3)**2)==F(1,19), 'axis tolerance')
    ell = F(9,4)
    require((2*ell-4)/(ell*ell+5)==F(8,161)<F(1,19), 'axis sign tolerance')
    require(root_parameters(11,F(13,6))==(F(1,37),F(1,35)), 'm=11 root parameters')
    require(root_parameters(12,F(761,348))==(F(4225,125329),F(7021,208275)), 'm=12 root parameters')
    require(F(7021,208275)<F(4225,125329), 'm=12 maximum tolerance')
    for n,p,q,expected in [(6,3,40,595),(7,1,19,680),(8,1,19,288),(10,1,24,495),
                          (14,1,35,1800),(13,4225,125329,1296974)]:
        require(slack(n,F(p,q))==expected>0, f'D({n};{p},{q})')
    for mu,last in [(F(3,40),6),(F(1,19),8),(F(1,24),10),(F(1,35),15),(F(4225,125329),13)]:
        require(slack(last,mu)>0 and slack(last+1,mu)<=0, 'strict admissibility endpoint')
    print('PASS: weights, tolerances and strict admissibility endpoints')
    require(F(156,256)==F(39,64), 'general eight-dimensional factor')
    check_table(args.article)
    print('All exact checks passed. No contact certificates or numerical searches were executed.')


if __name__ == '__main__':
    main()
