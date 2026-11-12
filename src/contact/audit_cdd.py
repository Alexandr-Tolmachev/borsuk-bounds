"""Second, slower replay of a random sample of leaves: the same cells, but every vertex set is enumerated by
cddlib in exact GMP arithmetic (no adjacency-closure shortcut, no pruning).
usage: audit_cdd.py cert.json.gz N [--workers W] [--seed S]"""
import os, sys, json, gzip, time, random
os.environ.setdefault('OMP_NUM_THREADS', '1')
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import verify as VF


def check_leaf_cdd(args):
    box, part, den = args
    import tri_cert as T
    from check_certificate import exact_vertices, max_sq_dist
    A = [[Fr(x, den) for x in row] for row in part[0]]; b = [Fr(x, den) for x in part[1]]
    balls, planes, m = T.outer(box)
    worst = 0.0; t0 = time.time()
    for i in range(6):
        tans, d, fok = T.build_cell(balls, planes, A, b, i)
        cell = [(a, bb, key[0]) for key, (a, bb) in tans.items()]
        if not all(T.tangent_ok(a, bb, balls[k]) for a, bb, k in cell): return box, False, None, 0
        H = list(planes) + T.cell_halfspaces(A, b, i) + [(a, bb) for a, bb, _ in cell]
        V = exact_vertices(H, 4)
        if V is None: return box, False, None, 0
        if not V: continue
        ex, fl = max_sq_dist(V)
        if ex is not None and ex >= 1: return box, False, float(ex) ** .5, 0
        worst = max(worst, float(fl) ** .5)
    return box, True, worst, time.time() - t0


if __name__ == '__main__':
    path, n = sys.argv[1], int(sys.argv[2])
    w = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 1
    seed = int(sys.argv[sys.argv.index('--seed') + 1]) if '--seed' in sys.argv else 2
    with gzip.open(path, 'rt') as f: data = json.load(f)
    L, vol = VF.leaves(data['tree'])
    jobs = [(box, data['partitions'][j], data['denominator']) for box, j in L if j is not None]
    jobs = random.Random(seed).sample(jobs, min(n, len(jobs)))
    from multiprocessing import Pool
    fails, worst, done = 0, 0.0, 0
    with Pool(w) as pool:
        for box, ok, wv, dt in pool.imap_unordered(check_leaf_cdd, jobs):
            done += 1
            if not ok: fails += 1; print('FAIL', box, wv, flush=True)
            else: worst = max(worst, wv)
            if done % 25 == 0: print(f'cdd audit {done}/{len(jobs)}, failures {fails}, worst {worst:.10f}', flush=True)
    print(json.dumps({'sampled': done, 'failures': fails, 'worst_diameter_bound': worst}), flush=True)
