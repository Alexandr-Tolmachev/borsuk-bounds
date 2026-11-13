# Novel upper bounds for Borsuk numbers

Certificates and verification programs for the article

> **Novel upper bounds for Borsuk numbers from contact points and spherical partitions**,
> arXiv preprint [arXiv:2611.16775](https://arxiv.org/abs/2611.16775).

The Borsuk number $b(n)$ is the least number of subsets of strictly smaller diameter needed to partition every bounded set of positive diameter in $\mathbb{R}^n$. The article proves the following upper bounds.

| $n$ | Upper bound for $b(n)$ | Proof |
|---:|---:|---|
| 4 | 6 | geometric certificate |
| 5 | 12 | geometric certificate |
| 6 | 31 | rational matrix certificate |
| 7 | 45 | analytic |
| 8 | 89 | analytic |
| 9 | 157 | analytic |
| 10 | 289 | analytic |
| 11 | 625 | analytic |
| 12 | 1245 | analytic |
| 13 | 2313 | analytic |
| 14 | 4977 | analytic |
| 15 | 7745 | analytic |
| 16 | 13729 | analytic |
| 17 | 24337 | analytic |
| 18 | 48673 | analytic |
| 19 | 97345 | analytic |

For $n \ge 9$ and nonnegative integers $a, b, c$ with $n-1 = 7a+8b+c$, the article proves

$$b(n) \le 1 + 88^a \cdot 156^b \cdot 2^c .$$

In particular, if $n-1 = 8q+r$ with $0 \le r < 8$, then $b(n) \le 1 + 156^q \cdot 2^r$. The geometric proofs of these analytic estimates are in the article. The arithmetic checker verifies their constants, algebraic identities and finite consequences.

## Quick verification

Use Python 3.10 or newer. Run all commands below from the repository root.

```sh
python3 verify.py inspect
python3 verify.py short
```

No packages need to be installed for these commands. They do not run the long geometric verifications.

- `inspect` checks the SHA-256 checksums of the code, data, certificates and logs (documentation files such as this README are not included), both subdivision trees, the rational covered volume, the power-diagram dimensions and the stored verification summaries. It also checks the outward correction of the printed diameter bounds. Reading a stored summary is not a geometric verification.
- `short` verifies the four rational matrix certificates for dimension six, including all 63 rays and 1007 bilinear products, and checks the algebraic identities, admissibility conditions and upper bounds for $6 \le n \le 19$ in exact arithmetic.

Both commands exit with status zero on success. Failed assertions or failed checks produce a nonzero exit status. Keep Python assertions enabled: do not use `-O`, `-OO` or `PYTHONOPTIMIZE`.

## Complete geometric verification

These commands rebuild and check every leaf of the certificates for dimensions four and five. They require NumPy, SciPy, and pycddlib 3 with its GMP backend.

On Debian or Ubuntu, install the build dependencies and create a Python environment:

```sh
sudo apt-get install build-essential python3-dev python3-venv libgmp-dev libcdd-dev
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -c "import numpy, scipy, cdd.gmp"
```

Then run:

```sh
python verify.py full --dimension 4 --workers 4
python verify.py full --dimension 5 --workers 4
```

The stored complete runs took about 75 minutes and 11.5 hours, respectively, using four processes. These times depend on the machine and numerical libraries. `--workers 1` uses one process.

A successful full run must end with the following values:

| Dimension | `checked` | `failures` | `covered_volume` | `complete_certificate` |
|---|---:|---:|---|---|
| 4 | 11819 | 0 | `"1"` | `true` |
| 5 | 6976 | 0 | `"1"` | `true` |

The checker regenerates tangent halfspaces, checks their support in rational arithmetic, certifies the complete vertex lists, and verifies strict diameter inequalities. Floating-point computations propose candidates or apply filters with proved error bounds; inconclusive cases use exact arithmetic. A failed numerical proposal is not a disproof of the bound: it means the verification did not finish successfully in that environment.

## Every verification command

| Command | What it checks | Dependencies |
|---|---|---|
| `python3 verify.py inspect` | file integrity, tree structure and stored summaries | standard library |
| `python3 verify.py six` | matrix and ray certificate for $b(6) \le 31$ | standard library |
| `python3 verify.py analytic` | constants, identities and analytic upper bounds | standard library |
| `python3 verify.py short` | `six` and `analytic` | standard library |
| `python verify.py full --dimension 4 --workers 4` | complete six-part geometric certificate | `requirements.txt` |
| `python verify.py full --dimension 5 --workers 4` | complete twelve-part geometric certificate | `requirements.txt` |
| `python verify.py cdd --dimension 4 --sample 120 --workers 4` | independent GMP vertex enumeration on 120 leaves | `requirements.txt` |
| `python verify.py cdd --dimension 5 --sample 120 --workers 4` | independent GMP vertex enumeration on 120 leaves | `requirements.txt` |

The independent `cdd` checks use a fixed random seed. They are additional sample checks, not replacements for complete verification. They can also take considerable time.

To run **all** stages, including both complete geometric checks and both independent samples:

```sh
python verify.py all --workers 4
```

This is a long computation. To save a run's output, redirect it to a new file, for example `python verify.py all --workers 4 > verification.log 2>&1`. The program never overwrites the stored certificate logs.

## Compare with the article

The repository is self-contained and does not require LaTeX. The article sources are distributed separately. If they are extracted into a directory `../article` containing `main.tex` and `sections/`, use

```sh
python3 verify.py short --article ../article
```

This additionally checks that the matrices and rays in the appendix agree with `data/six.json`, and that the numerical bounds in Table 1 agree with the arithmetic checker. The path may be absolute.

## Files

| Path | Contents |
|---|---|
| `verify.py` | entry point for all verification commands |
| `checks/` | exact matrix, arithmetic and structural checks |
| `src/contact/` | rational enclosures, tangent construction and geometric verifiers |
| `certificates/b4_k6/` | four-dimensional tree and complete/sample verification logs |
| `certificates/b5_k12/` | five-dimensional tree and complete/sample verification logs |
| `data/six.json` | matrices and integer ray generators for the dimension-six proof |
| `requirements.txt` | dependencies for geometric verification |
| `SHA256SUMS` | SHA-256 checksums of the code, data, certificates and logs |
| `THIRD_PARTY_NOTICES.md` | attribution for external algorithms and code |

Each `cert.json.gz` is gzip-compressed JSON. The field `tree` is a preorder binary tree: `0,1,2,3` split that parameter coordinate at its midpoint; `-(j+1)` denotes a leaf using diagram `j`. Each diagram consists of integer site vectors and offsets over the common denominator $10^9$. Both trees cover the full parameter cube $[0,1]^4$. The rational outer bodies and tangent halfspaces are reconstructed during verification.

The files `verify_full.log` and `audit_cdd.log` record complete and independent sample runs, respectively. Their role and the error bounds used by the verifier are explained in the article's verification appendix.
