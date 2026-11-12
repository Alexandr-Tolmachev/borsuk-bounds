"""Complete geometric verification in dimension four.

Rebuilds each dyadic box and its rational outer model, proposes tangent
halfspaces, verifies their support and checks all cell diameters.
Run through the repository-level verify.py entry point.
"""
import os, sys, json, gzip, time, random
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))

FORMAT = 'b4-six-parts-v1'
ROOT = ((0, 0),) * 4




def leaves(tree):
    """parse the preorder tree; returns [(box, partition index)] and the exact covered volume"""
    out, stack, pos, vol = [], [ROOT], 0, Fr(0)
    while stack:
        assert pos < len(tree), 'truncated tree'
        box = stack.pop(); t = tree[pos]; pos += 1
        if t is None:
            out.append((box, None)); continue
        assert type(t) is int, 'bad token'
        if t >= 0:
            assert t < 4
            k, l = box[t]; L, R = list(box), list(box); L[t] = (2 * k, l + 1); R[t] = (2 * k + 1, l + 1)
            stack += [tuple(R), tuple(L)]
        else:
            out.append((box, -t - 1)); vol += Fr(1, 2 ** sum(l for _, l in box))
    assert pos == len(tree), 'trailing tokens'
    return out, vol


def check_leaf(args):
    box, part, den = args
    import tri_cert as T
    A = [[Fr(x, den) for x in row] for row in part[0]]; b = [Fr(x, den) for x in part[1]]
    assert len(A) == 6 and all(len(r) == 4 for r in A) and len(b) == 6
    t0 = time.time()
    balls, planes, m = T.outer(box)
    worst = 0.0
    for i in range(6):
        tans, d, fok = T.build_cell(balls, planes, A, b, i)
        cell = [(a, bb, key[0]) for key, (a, bb) in tans.items()]
        ok, w, nv = T.check_cell(balls, planes, A, b, i, cell)
        if not ok: return box, False, w, time.time() - t0
        worst = max(worst, w)
    return box, True, worst, time.time() - t0


def verify(path, workers=1, sample=None):
    with gzip.open(path, 'rt') as f: data = json.load(f)
    assert data['format'] == FORMAT and data['r_interval'] == ['1/4', '2/5']
    L, vol = leaves(data['tree'])
    unresolved = [x for x in L if x[1] is None]
    jobs = [(box, data['partitions'][j], data['denominator']) for box, j in L if j is not None]
    if sample: jobs = random.Random(1).sample(jobs, min(sample, len(jobs)))
    print(f'{len(L)} leaves, {len(unresolved)} unresolved, covered volume {float(vol):.9f}, checking {len(jobs)}', flush=True)
    fails, worst, done, t0 = [], 0.0, 0, time.time()
    from multiprocessing import Pool
    with Pool(workers) as pool:
        for box, ok, w, dt in pool.imap_unordered(check_leaf, jobs, chunksize=4):
            done += 1
            if not ok: fails.append(box); print('FAIL', box, w, flush=True)
            else: worst = max(worst, w)
            if done % 500 == 0: print(f'checked {done}/{len(jobs)}, worst {worst:.6f} [{time.time()-t0:.0f}s]', flush=True)
    complete = not fails and not unresolved and vol == 1 and not sample
    print(json.dumps({'checked': done, 'failures': len(fails), 'worst_diameter_bound': worst,
                      'covered_volume': str(vol), 'complete_certificate': complete}), flush=True)
    return complete


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate')
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    if not __debug__:
        parser.error('Python assertions must be enabled.')
    if args.workers < 1:
        parser.error('--workers must be positive')
    sys.exit(0 if verify(args.certificate, args.workers) else 1)
