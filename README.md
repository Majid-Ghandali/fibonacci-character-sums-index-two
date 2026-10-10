# Fibonacci Character Sums at Rank Index Two

[![Manuscript DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23264457.svg)](https://doi.org/10.5281/zenodo.23264457)
[![Code Archive DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23263042.svg)](https://doi.org/10.5281/zenodo.23263042)
[![License: MIT (code)](https://img.shields.io/badge/License-MIT%20(code)-yellow.svg)](LICENSE)
[![License: CC BY 4.0 (paper)](https://img.shields.io/badge/License-CC%20BY%204.0%20(paper)-lightgrey.svg)](Paper/LICENSE)
[![Python](https://img.shields.io/badge/Python-3.12.7-blue.svg)](https://www.python.org/)

**Majid Ghandali** Â· Independent Researcher, Tehran, Iran Â· 2026<br>
ORCID: [0009-0001-1097-1770](https://orcid.org/0009-0001-1097-1770)

Reproducibility materials for the preprint:

> *Fibonacci Character Sums at Rank Index Two and Jacobsthal Sums*

> **Proof versus computation.** The mathematical results are proved in the manuscript. Computations in this repository provide finite-range verification of the stated identities and reported numerical results; they are not a substitute for the proofs.

---

## Citation and DOI records

The manuscript and the computational release are separate research objects:

| Object | Use | DOI |
| --- | --- | --- |
| Manuscript preprint | Mathematical results and proofs | [10.5281/zenodo.23264457](https://doi.org/10.5281/zenodo.23264457) |
| Reproducibility archive v1.0.1 | Software, data, and reproducibility reruns | [10.5281/zenodo.23263042](https://doi.org/10.5281/zenodo.23263042) |

- Manuscript citation metadata: [`Paper/CITATION.cff`](Paper/CITATION.cff)
- Software citation metadata: [`CITATION.cff`](CITATION.cff)

---

## Quick start

```bash
git clone https://github.com/Majid-Ghandali/fibonacci-character-sums-index-two.git
cd fibonacci-character-sums-index-two

python -m venv .venv
source .venv/bin/activate          # macOS/Linux
# .\.venv\Scripts\Activate.ps1     # Windows PowerShell

python -m pip install --upgrade pip
python -m pip install -r Requirements.txt
```

Run a small smoke verification:

```bash
python code/ghandali_verify_index2.py \
  --limit 10000 \
  --sum-limit 10000 \
  --jac-limit 10000 \
  --workers 1 \
  --sum-indices 2,4 \
  --out results/reproduced-smoke-N10000 \
  --fresh
```

A successful run ends with:

```text
VERDICT: PASS
theorem mismatches: 0
```

Windows PowerShell:

```powershell
python code\ghandali_verify_index2.py --limit 10000 --sum-limit 10000 --jac-limit 10000 --workers 1 --sum-indices 2,4 --out results\reproduced-smoke-N10000 --fresh
```

---

## Main results (statements only)

For a prime $p\ne 2,5$, let

$$
\chi_p=\left(\frac{\cdot}{p}\right),\qquad
\varepsilon_p=\left(\frac{5}{p}\right),\qquad
n_p=p-\varepsilon_p,
$$

and let $z(p)$ be the rank of apparition of $p$ in the Fibonacci sequence. The rank index is

$$
i_p=\frac{n_p}{z(p)}.
$$

The Fibonacci character sum studied in the paper is

$$
S(p)=\sum_{j=1}^{z(p)}\chi_p(F_j).
$$

The paper focuses on the rank-index-two case

$$
i_p=2.
$$

### Split case

If $\left(\frac5p\right)=1$ and $i_p=2$, then $p\equiv1\pmod4$. Write


$$p=a^2+b^2,\qquad a\ \text{odd},\qquad b\ \text{even},\qquad b>0,$$


and let $$a^*\in\{\pm a\}$$ be chosen by

$$
a^*\equiv1\pmod4.
$$

Let

$$
\delta_p=
\begin{cases}
1,&p\equiv1\pmod5,\\
-1,&p\equiv4\pmod5.
\end{cases}
$$

Then

$$
S(p)=
\begin{cases}
\delta_p a^*,&p\equiv1\pmod8,\\
-\delta_p,&p\equiv5\pmod8.
\end{cases}
$$

### Inert case

If $\left(\frac5p\right)=-1$ and $i_p=2$, then $p\equiv1\pmod4$, and

$$
S(p)=
\begin{cases}
\frac12\mathcal J(5)\in\{\pm b\},&p\equiv1\pmod8,\\
0,&p\equiv5\pmod8,
\end{cases}
$$

where

$$
\mathcal J(c)=\sum_{x\in\mathbb F_p}\chi_p(x)\chi_p(x^2+c)
$$

is the classical Jacobsthal sum.

For the full hypotheses, reductions, lemmas, and proofs, see [`Paper/ghandali_index2.tex`](Paper/ghandali_index2.tex).

---

## Verification headline (committed evidence)

Reported release run through $p\le 2{,}000{,}000$:

| Quantity | Value |
| --- | ---: |
| Prime limit $N$ | 2,000,000 |
| Primes in rank inventory | 148,931 |
| Direct $S(p)$ and Jacobsthal record rows selected | 52,477 |
| Split index-two primes, $p\equiv1\pmod8$ | 6,551 |
| Split index-two primes, $p\equiv5\pmod8$ | 6,583 |
| Inert index-two primes, $p\equiv1\pmod8$ | 14,591 |
| Inert index-two primes, $p\equiv5\pmod8$ | 14,696 |
| Theorem-level mismatches | **0** |
| Exploratory mismatches | **0** |

Authoritative machine-readable summary:

```text
results/release-v1.0.1/verification-summary.json
```

Run log:

```text
results/release-v1.0.1/verification-run.log
```

These are finite computational results offered as consistency evidence, not as numerical proofs.

---

## Verification scope

The verifier separates three layers:

| Layer | Role |
| --- | --- |
| Theorem-level checks | Finite-range verification of identities proved in the manuscript |
| Exploratory checks | Open-pattern checks and auxiliary observations not claimed as theorems |
| Density context | Numerical ratios recorded as context only, not finite proofs of density statements |

The density computations in the output files are not used as proofs. In particular, the inert density statement in the manuscript is recorded as a GRH-conditional corollary from the work of Roskam and Chen, not as a conclusion drawn from finite computation.

---

## What the program computes

| Symbol / quantity | Meaning |
| --- | --- |
| $z(p)$ | Rank of apparition of $p$ in the Fibonacci sequence |
| $i_p$ | Rank index $(p-(5/p))/z(p)$ |
| $S(p)$ | Fibonacci character sum $\sum_{j=1}^{z(p)}\chi_p(F_j)$ |
| $\mathcal J(c)$ | Jacobsthal sum $\sum_{x\in\mathbb F_p}\chi_p(x)\chi_p(x^2+c)$ |
| $a^*,b$ | Two-square parameters for $p=a^2+b^2$, $a^*\equiv1\pmod4$ |
| $\delta_p$ | Split sign determined by $p\bmod 5$ |
| density ratios | Contextual numerical ratios relative to Artin's constant |

The program records rank-index inventories, direct Fibonacci character sums for selected rank indices, direct Jacobsthal sums where required, theorem-level mismatch counts, exploratory mismatch counts, numerical tables, and contextual density data.

---

## Full verification command

The archived release run used the public release script:

```bash
python code/ghandali_verify_index2.py \
  --limit 2000000 \
  --sum-limit 2000000 \
  --jac-limit 2000000 \
  --workers 14 \
  --sum-indices 2,4 \
  --out results/release-v1.0.1 \
  --fresh
```

Windows PowerShell:

```powershell
python code\ghandali_verify_index2.py --limit 2000000 --sum-limit 2000000 --jac-limit 2000000 --workers 14 --sum-indices 2,4 --out results\release-v1.0.1 --fresh
```

The release run completed with:

```text
VERDICT: PASS
theorem mismatches: 0
exploratory mismatches: 0
```

For new local reruns, always use a new output directory, for example:

```bash
python code/ghandali_verify_index2.py \
  --limit 2000000 \
  --sum-limit 2000000 \
  --jac-limit 2000000 \
  --workers 14 \
  --sum-indices 2,4 \
  --out results/reproduced-N2000000 \
  --fresh
```

Do **not** overwrite `results/release-v1.0.1/`.

Full CLI help:

```bash
python code/ghandali_verify_index2.py --help
```

---

## Frozen reproducibility release (v1.0.1)

Canonical software release: **v1.0.1**<br>
Executable artifact: [`code/ghandali_verify_index2.py`](code/ghandali_verify_index2.py)

Official evidence bundle:

```text
results/release-v1.0.1/
```

This directory is checksum-recorded release evidence. Do **not** rename files in it, normalize its line endings, edit it, or overwrite it with a new rerun.

Authoritative records:

```text
results/release-v1.0.1/SHA256SUMS.txt
results/release-v1.0.1/verification-summary.json
results/release-v1.0.1/verification-run.log
```

Additional archived outputs:

```text
results/release-v1.0.1/rank-index-all-primes.csv
results/release-v1.0.1/ghandali-index2-records.csv
results/release-v1.0.1/numerical-table.csv
results/release-v1.0.1/numerical-table.tex
results/release-v1.0.1/verification-workbook.xlsx
```

Fresh local smoke test:

```bash
python code/ghandali_verify_index2.py \
  --limit 10000 \
  --sum-limit 10000 \
  --jac-limit 10000 \
  --workers 1 \
  --sum-indices 2,4 \
  --out results/reproduced-smoke-N10000 \
  --fresh
```

### Integrity verification

From the repository root, if GNU `sha256sum` is available:

```bash
sha256sum -c results/release-v1.0.1/SHA256SUMS.txt
```

Windows PowerShell example for a single file:

```powershell
Get-FileHash code\ghandali_verify_index2.py -Algorithm SHA256
# Compare the digest with the corresponding line in results\release-v1.0.1\SHA256SUMS.txt
```

The official executable digest and output-file digests are recorded in the checksum manifest. Do not clean or rewrite release evidence to satisfy cosmetic checks; integrity is defined by `SHA256SUMS.txt`.

---

## Repository layout

```text
Paper/                               manuscript files (CC BY 4.0): TeX, BibTeX, PDF, LICENSE, CITATION.cff
code/ghandali_verify_index2.py        frozen executable for release v1.0.1
results/release-v1.0.1/               frozen, checksum-recorded release evidence
results/reproduced-*/                 optional local reruns, not part of the frozen release
docs/reproducibility.md               optional reproduction notes
LICENSE                               MIT license for code and computational materials
Paper/LICENSE                         CC BY 4.0 license for manuscript text
LICENSES.md                           scope of the two licenses
CITATION.cff                          software citation metadata
Requirements.txt                      Python environment notes
README.md                             repository overview
```

---

## Requirements

- Tested release environment: **Python 3.12.7**
- The verifier is written in Python.
- If `Requirements.txt` lists no third-party packages, the verifier uses only the Python standard library.

Install the repository requirements, if any:

```bash
python -m pip install -r Requirements.txt
```

---

## Manuscript compilation

Source: [`Paper/ghandali_index2.tex`](Paper/ghandali_index2.tex)

```bash
cd Paper
pdflatex -interaction=nonstopmode ghandali_index2.tex
bibtex ghandali_index2
pdflatex -interaction=nonstopmode ghandali_index2.tex
pdflatex -interaction=nonstopmode ghandali_index2.tex
```

A TeX Live or MiKTeX installation providing the packages declared in the source is required.

---

## BibTeX

**Manuscript**

```bibtex
@misc{Ghandali2026IndexTwoPreprint,
  author       = {Ghandali, Majid},
  title        = {Fibonacci Character Sums at Rank Index Two and Jacobsthal Sums},
  year         = {2026},
  version      = {v1.0.0-preprint},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23264457},
  url          = {https://doi.org/10.5281/zenodo.23264457}
}
```

**Software / reproducibility archive**

```bibtex
@misc{Ghandali2026IndexTwoArchive,
  author       = {Ghandali, Majid},
  title        = {Reproducibility Archive for Fibonacci Character Sums at Rank Index Two and Jacobsthal Sums},
  year         = {2026},
  version      = {1.0.1},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23263042},
  url          = {https://doi.org/10.5281/zenodo.23263042}
}
```

---

## License

This repository uses a deliberate dual-license structure:

| Content | Path | License |
| --- | --- | --- |
| Code, scripts, verification outputs, and computational artifacts | Root [`LICENSE`](LICENSE); `code/`, `results/`, â€¦ | **MIT** |
| Manuscript text, bibliography, and PDF | [`Paper/`](Paper/); [`Paper/LICENSE`](Paper/LICENSE) | **CC BY 4.0** |

Full scope statement: [`LICENSES.md`](LICENSES.md).

---

## Author

**Majid Ghandali**<br>
Independent Researcher, Tehran, Iran<br>
Email: [majid.ghandali@gmail.com](mailto:majid.ghandali@gmail.com)<br>
ORCID: [0009-0001-1097-1770](https://orcid.org/0009-0001-1097-1770)
