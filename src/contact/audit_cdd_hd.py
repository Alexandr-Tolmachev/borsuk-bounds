"""Second, slower replay of a random sample of leaves of a k-part certificate in R^n: the same cells, but every
vertex set is enumerated by cddlib in exact GMP arithmetic (no adjacency-closure shortcut, no pruning).
usage: audit_cdd_hd.py cert.json.gz N [--workers W] [--seed S] [--done audit_done.jsonl]"""
import os, sys, json, gzip, time, random
os.environ.setdefault('OMP_NUM_THREADS', '1')
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import verify_hd as VF


def check_leaf_cdd(args):
    box, part, den, n, k, symmetry = args
    import hd_cert as T
    from check_certificate import exact_vertices, max_sq_dist
    T.set_dim(n); T.MIN_SEP = 2e-3
    A = [[Fr(x, den) for x in row] for row in part[0]]; b = [Fr(x, den) for x in part[1]]
    cells = range(k)
    if symmetry is not None:
        import sym as SY
        SIG = [[int(x) for x in s_] for s_ in SY.signs(symmetry, n)]; g = len(SIG)
        for j in range(0, k, g):
            for t in range(g):
                assert A[j + t] == [SIG[t][c] * A[j][c] for c in range(n)] and b[j + t] == b[j]
        cells = range(0, k, g)
    balls, planes, m = T.outer(box)
    worst = 0.0; t0 = time.time()
    def one(i):
        tans, d, fok = T.build_cell(balls, planes, A, b, i)
        cell = [(a, bb, key[0]) for key, (a, bb) in tans.items()]
        if not all(T.tangent_ok(a, bb, balls[kk]) for a, bb, kk in cell): return False, None
        H = list(planes) + T.cell_halfspaces(A, b, i) + [(a, bb) for a, bb, _ in cell]
        V = exact_vertices(H, n)
        if V is None: return False, None
        if not V: return True, 0.0
        ex, fl = max_sq_dist(V)
        if ex is not None and ex >= 1: return False, float(ex) ** .5
        return True, float(fl) ** .5
    for i in cells:
        for sep in (2e-3, 0.0):          # the two tangent constructions of verify_hd
            T.MIN_SEP = sep
            ok, w = one(i)
            if ok: break
        if not ok: return box, False, w, 0
        worst = max(worst, w)
    return box, True, worst, time.time() - t0


if __name__ == '__main__':
    path, N = sys.argv[1], int(sys.argv[2])
    w = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 1
    seed = int(sys.argv[sys.argv.index('--seed') + 1]) if '--seed' in sys.argv else 2
    with gzip.open(path, 'rt') as f: data = json.load(f)
    n, k = data['dimension'], data['parts']
    L, vol = VF.leaves(data['tree'])
    jobs = [(box, data['partitions'][j], data['denominator'], n, k, data.get('symmetry')) for box, j in L if j is not None]
    jobs = random.Random(seed).sample(jobs, min(N, len(jobs)))
    # optional --done FILE: per-leaf results are appended there and leaves already in it are not redone
    # (lets an interrupted audit resume on the same sample)
    done_file = sys.argv[sys.argv.index('--done') + 1] if '--done' in sys.argv else None
    prev = {}
    if done_file and os.path.exists(done_file):
        for line in open(done_file):
            d = json.loads(line); prev[json.dumps(d['box'])] = d
    key = lambda box: json.dumps([list(x) for x in box])
    from multiprocessing import Pool
    fails, worst, done = 0, 0.0, 0
    for d in prev.values():
        if key(d['box']) not in {key(j[0]) for j in jobs}: continue
        done += 1
        if not d['ok']: fails += 1
        else: worst = max(worst, d['worst'])
    todo = [j for j in jobs if key(j[0]) not in prev]
    if prev: print(f'resuming: {done} leaves already audited, {len(todo)} to go', flush=True)
    with Pool(w) as pool:
        for box, ok, wv, dt in pool.imap_unordered(check_leaf_cdd, todo):
            done += 1
            if done_file:
                with open(done_file, 'a') as f: f.write(json.dumps({'box': [list(x) for x in box], 'ok': ok, 'worst': wv}) + '\n')
            if not ok: fails += 1; print('FAIL', box, wv, flush=True)
            else: worst = max(worst, wv)
            if done % 10 == 0: print(f'cdd audit {done}/{len(jobs)}, failures {fails}, worst {worst:.10f}', flush=True)
    print(json.dumps({'sampled': done, 'failures': fails, 'worst_diameter_bound': worst}), flush=True)
