"""Check file integrity, certificate structure and stored verification summaries."""
from fractions import Fraction
from pathlib import Path
import gzip
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def inspect_certificate(n, k, expected_leaves, expected_diagrams):
    folder = ROOT / 'certificates' / f'b{n}_k{k}'
    with gzip.open(folder / 'cert.json.gz', 'rt') as stream:
        data = json.load(stream)
    require(data['denominator'] == 10**9, 'unexpected denominator')
    require(data['format'] == ('b4-six-parts-v1' if n == 4 else 'three-contact-partition-v1'), 'format')
    if n == 5:
        require(data['dimension'] == n and data['parts'] == k, 'dimension/parts')
        require(data.get('symmetry') is None, 'this audit expects the unsymmetrized certificate')
    require(data['r_interval'] == (['1/4', '2/5'] if n == 4 else ['1/4', '5/12']), 'radius interval')
    parts = data['partitions']
    require(len(parts) == expected_diagrams, 'diagram count')
    for a, b in parts:
        require(len(a) == k and len(b) == k and all(len(row) == n for row in a), 'diagram shape')
        require(all(type(t) is int for row in a for t in row) and all(type(t) is int for t in b), 'noninteger data')
    stack, leaves, volume, depths = [(0, 0, 0, 0)], 0, Fraction(0), []
    for token in data['tree']:
        require(bool(stack), 'trailing tree tokens')
        depth = stack.pop()
        require(type(token) is int, 'unresolved or invalid tree token')
        if token >= 0:
            require(token < 4, 'invalid split coordinate')
            child = tuple(t + (j == token) for j, t in enumerate(depth))
            stack.extend((child, child))
        else:
            require(0 <= -token-1 < len(parts), 'invalid diagram index')
            leaves += 1
            depths.append(sum(depth))
            volume += Fraction(1, 2**sum(depth))
    require(not stack and leaves == expected_leaves and volume == 1, 'incomplete subdivision')
    records = [json.loads(line) for line in (folder / 'verify_full.log').read_text().splitlines() if line.startswith('{')]
    report = records[-1]
    require(report['checked'] == leaves and report['failures'] == 0 and report['covered_volume'] == '1'
            and report['complete_certificate'] is True, 'stored full log does not match certificate')
    from fractions import Fraction as F
    exact_report = json.loads([line for line in (folder/'verify_full.log').read_text().splitlines()
                               if line.startswith('{')][-1], parse_float=F)
    target = F('0.9999999994') if n==4 else F('0.99999995')
    require(exact_report['worst_diameter_bound']**2+F(2,10**12)<target**2,
            'outward diameter correction exceeds the published bound')
    return {'dimension': n, 'parts': k, 'leaves': leaves, 'diagrams': len(parts),
            'depth': [min(depths), max(depths)], 'volume': str(volume),
            'stored_full_verification': report, 'geometric_replay_performed': False}


def main():
    verified = 0
    for line in (ROOT / 'SHA256SUMS').read_text().splitlines():
        digest, relative = line.split('  ', 1)
        path = ROOT / relative
        require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, f'checksum mismatch: {relative}')
        verified += 1
    results = [inspect_certificate(4, 6, 11819, 10615), inspect_certificate(5, 12, 6976, 5954)]
    print(json.dumps({'files_verified': verified, 'certificates': results,
                      'scope': 'checksums, structure, stored logs; no geometric replay'}, indent=2))


if __name__ == '__main__':
    main()
