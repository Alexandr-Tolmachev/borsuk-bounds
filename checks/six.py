"""Exact verification of the dimension-six matrix certificate."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd, lcm
from pathlib import Path
import re
import json
import argparse
import sys
from time import perf_counter


def require(condition, message):
    if not condition:
        raise ValueError(message)


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def determinant(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    if n == 0:
        return 1
    sign, previous = 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = a[i][j]*a[k][k] - a[i][k]*a[k][j]
                require(value % previous == 0, 'nonexact Bareiss division')
                a[i][j] = value // previous
        previous = a[k][k]
    return sign*a[-1][-1]


def primitive(row):
    g = reduce(gcd, map(abs, row))
    require(g > 0, 'zero row')
    return tuple(x//g for x in row)


def extreme_rays(a):
    require(any(determinant([a[i] for i in indices]) != 0
                for indices in combinations(range(len(a)), 5)), 'cone not pointed')
    result = set()
    for indices in combinations(range(len(a)), 4):
        matrix = [a[i] for i in indices]
        row = [(-1)**k * determinant([r[:k] + r[k+1:] for r in matrix])
               for k in range(5)]
        if not any(row):
            continue
        row = primitive(row)
        if all(dot(r, row) >= 0 for r in a):
            result.add(row)
        if all(dot(r, row) <= 0 for r in a):
            result.add(tuple(-x for x in row))
    return sorted(result)


# All functional vectors are scaled by 20, which leaves their normal cones unchanged.
def configuration():
    eye = [[int(i == j) for j in range(5)] for i in range(5)]
    vertices = [tuple(44*x for x in e) for e in eye]
    vertices += [tuple(-40*x for x in e) for e in eye]
    for signs in product((-1, 1), repeat=5):
        negatives = signs.count(-1)
        if negatives in (0, 1, 2, 4):
            weight = 19 if negatives == 4 else 20
            vertices.append(tuple(weight*x for x in signs))
    require(len(vertices) == len(set(vertices)) == 31, 'incorrect configuration')
    return vertices



def verify_case(case, vertices):
    name, v, groups, denominator, numerator, a, expected_count, expected_min = case
    require(v in vertices, 'missing representative')
    allowed = {primitive([x-y for x, y in zip(v, w)]) for w in vertices if w != v}
    permutation = list(range(5))
    for group in groups:
        for i, j in zip(group[:-1], group[1:]):
            row = [0]*5
            row[i], row[j] = -1, 1
            allowed.add(tuple(row))
        for i, j in zip(group, reversed(group)):
            permutation[i] = j
    require(all(primitive(row) in allowed for row in a), 'unjustified cone inequality')
    e = [[F(x, denominator) for x in row] for row in numerator]
    s = [[int(i == j) + e[i][j] for j in range(5)] for i in range(5)]
    require(all(e[i][j] == e[j][i] for i in range(5) for j in range(5)), 'not symmetric')
    # Sylvester's criterion: ||E||_op < 3/40, with exact integer determinants.
    for sign in (-1, 1):
        matrix = [[3*denominator*int(i == j) + sign*40*numerator[i][j]
                   for j in range(5)] for i in range(5)]
        require(all(determinant([row[:k] for row in matrix[:k]]) > 0
                    for k in range(1, 6)), 'operator norm bound failed')
    # Invariance and a positive coefficient of the ordinary scalar product
    # within each coordinate block justify the rearrangement argument.
    for group in groups:
        diagonal = {s[i][i] for i in group}
        off = {s[i][j] for i in group for j in group if i != j}
        require(len(diagonal) == 1 and len(off) <= 1, 'invalid symmetric block')
        if off:
            require(next(iter(diagonal)) > next(iter(off)), 'negative rearrangement coefficient')
    for g in groups:
        for h in groups:
            if g != h:
                require(len({s[i][j] for i in g for j in h}) == 1, 'invalid cross block')
    rays = extreme_rays(a)
    require(len(rays) == expected_count, 'incorrect ray count')
    values = [sum(F(r[i])*s[i][j]*t[permutation[j]] for i in range(5) for j in range(5))
              for r in rays for t in rays]
    require(min(values) == expected_min and min(values) >= 0, 'bilinear certificate failed')
    print(f'{name}: {len(a)} inequalities, {len(rays)} rays, '
          f'{len(values)} products, minimum={min(values)}; PASS')
    return rays


def load_data(paths):
    matrices = {}
    pattern = re.compile(r"([AER])_\{([+\-24])\}=(?:\\frac1\{(\d+)\})?\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", re.S)
    for path in dict.fromkeys(paths):
        for name, typ, denominator, body in pattern.findall(path.read_text()):
            rows = [[F(int(x),int(denominator or 1)) for x in row.split('&')]
                    for row in body.strip().split('\\\\') if row.strip()]
            require(all(len(row)==5 for row in rows), 'matrix width in LaTeX')
            require((name,typ) not in matrices, 'duplicate matrix in LaTeX')
            matrices[name,typ] = rows
    require(set(matrices)=={(name,typ) for name in 'AER' for typ in '+-24'}, 'missing or extra matrix')
    descriptions = [('+','positive_axis',(44,0,0,0,0),[[0],[1,2,3,4]],17,F(89,100)),
                    ('-','negative_axis',(-40,0,0,0,0),[[0],[1,2,3,4]],15,F(53,50)),
                    ('2','two_negative_signs',(-20,-20,20,20,20),[[0,1],[2,3,4]],18,F(12,125)),
                    ('4','four_negative_signs',(-19,-19,-19,-19,19),[[0,1,2,3],[4]],13,F(209,100))]
    for typ,name,v,groups,count,minimum in descriptions:
        a,e,r = (matrices[key,typ] for key in 'AER')
        require(len(e)==5 and all(x.denominator==1 for rows in (a,r) for row in rows for x in row),
                'invalid matrix data')
        denominator = lcm(*(x.denominator for row in e for x in row))
        numerator = [[int(x*denominator) for x in row] for row in e]
        rays = [tuple(int(x) for x in row) for row in r]
        require(len(rays)==len(set(rays))==count, 'duplicate or missing printed ray')
        case = (name,v,groups,denominator,numerator,[[int(x) for x in row] for row in a],count,minimum)
        yield case,set(rays)


def main():
    start = perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--article', type=Path, help='also compare with the separate LaTeX source directory')
    args = parser.parse_args()
    path = Path(__file__).resolve().parents[1]/'data/six.json'
    data = json.loads(path.read_text())
    cases = []
    for row in data:
        rays = [tuple(x) for x in row['rays']]
        require(len(rays)==len(set(rays))==row['ray_count'], 'duplicate or missing certificate ray')
        case = (row['name'],tuple(row['representative']),row['groups'],row['E_denominator'],
                row['E_numerator'],row['A'],row['ray_count'],F(row['bilinear_minimum']))
        cases.append((case,set(rays)))
    require([c[0][0] for c in cases]==['positive_axis','negative_axis','two_negative_signs','four_negative_signs'],
            'missing or repeated region class')
    if args.article:
        printed = list(load_data([args.article/'sections/rays.tex']))
        require(cases==printed, 'certificate data differ from the LaTeX appendix')
        print('PASS: certificate data agree with the LaTeX appendix')
    vertices = configuration()
    row_subsets = pair_count = 0
    for case,printed_rays in cases:
        rays = verify_case(case,vertices)
        require(set(rays)==printed_rays, 'computed rays differ from the article')
        from math import comb
        row_subsets += comb(len(case[5]),4)
        pair_count += len(rays)**2
    require(row_subsets==301 and pair_count==1007, 'printed verification counts')
    print(f'PASS: dimension-six certificate, {row_subsets} active row subsets, {pair_count} ray products')
    print(f'All exact matrix checks passed in {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
