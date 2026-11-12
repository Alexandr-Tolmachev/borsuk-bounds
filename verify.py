#!/usr/bin/env python3
"""Verify the certificates and arithmetic for upper bounds on Borsuk numbers."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run(path, *args):
    subprocess.run([sys.executable, str(ROOT/path), *map(str, args)], cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['inspect', 'short', 'six', 'analytic', 'full', 'cdd', 'all'],
                        help='all includes both full geometric verifications and the independent cdd samples')
    parser.add_argument('--dimension', type=int, choices=[4, 5], help='required for full and cdd')
    parser.add_argument('--workers', type=int, default=1, help='parallel processes for full and cdd (default: 1)')
    parser.add_argument('--sample', type=int, default=120, help='leaves per independent cdd sample (default: 120)')
    parser.add_argument('--article', type=Path, help='optional LaTeX source directory for short, six, analytic or all')
    args = parser.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE', '') not in ('', '0'):
        parser.error('Assertions must be enabled: remove -O, -OO and PYTHONOPTIMIZE.')
    if args.workers < 1 or args.sample < 1:
        parser.error('--workers and --sample must be positive')
    if args.mode in ('full', 'cdd') and args.dimension is None:
        parser.error('--dimension 4 or 5 is required for full and cdd')
    if args.dimension is not None and args.mode not in ('full', 'cdd'):
        parser.error('--dimension is only used with full and cdd')
    if args.article and args.mode not in ('short', 'six', 'analytic', 'all'):
        parser.error('--article is only used with short, six, analytic or all')
    if args.article and not all((args.article/name).is_file()
                                for name in ('sections/rays.tex', 'sections/introduction.tex')):
        parser.error(f'{args.article} does not contain the article sources (sections/rays.tex, sections/introduction.tex)')
    article_args = ['--article', args.article.resolve()] if args.article else []

    def short():
        run('checks/six.py', *article_args)
        run('checks/analytic.py', *article_args)

    def geometric(n):
        script = 'verify.py' if n == 4 else 'verify_hd.py'
        certificate = ROOT/'certificates'/('b4_k6' if n == 4 else 'b5_k12')/'cert.json.gz'
        run(Path('src/contact')/script, certificate, '--workers', args.workers)

    def cdd(n):
        script = 'audit_cdd.py' if n == 4 else 'audit_cdd_hd.py'
        certificate = ROOT/'certificates'/('b4_k6' if n == 4 else 'b5_k12')/'cert.json.gz'
        command = [sys.executable, str(ROOT/'src/contact'/script), str(certificate), str(args.sample),
                   '--workers', str(args.workers)]
        report = None
        with subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE, text=True) as process:
            for line in process.stdout:
                print(line, end='', flush=True)
                if line.startswith('{'):
                    report = json.loads(line)
            code = process.wait()
        count = min(args.sample, 11819 if n == 4 else 6976)
        if code or report is None or report.get('failures') != 0 or report.get('sampled') != count:
            raise SystemExit('Independent cdd verification failed or did not check the requested sample.')

    if args.mode == 'inspect':
        run('checks/structure.py')
    elif args.mode == 'short':
        short()
    elif args.mode in ('six', 'analytic'):
        run(f'checks/{args.mode}.py', *article_args)
    elif args.mode == 'full':
        run('checks/structure.py')
        geometric(args.dimension)
    elif args.mode == 'cdd':
        run('checks/structure.py')
        cdd(args.dimension)
    else:
        run('checks/structure.py')
        short()
        geometric(4)
        geometric(5)
        cdd(4)
        cdd(5)


if __name__ == '__main__':
    main()
