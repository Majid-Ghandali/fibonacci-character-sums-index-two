# Fibonacci Character Sums at Rank Index Two

[![Manuscript DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23264456.svg)](https://doi.org/10.5281/zenodo.23264456)
[![Code Archive DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23263041.svg)](https://doi.org/10.5281/zenodo.23263041)
[![License: MIT (code)](https://img.shields.io/badge/License-MIT%20(code)-yellow.svg)](LICENSE)
[![License: CC BY 4.0 (paper)](https://img.shields.io/badge/License-CC%20BY%204.0%20(paper)-lightgrey.svg)](Paper/LICENSE)
[![Python](https://img.shields.io/badge/Python-3.12.7-blue.svg)](https://www.python.org/)

**Majid Ghandali** Independent Researcher, Tehran, Iran 2026.

ORCID: [0009-0001-1097-1770](https://orcid.org/0009-0001-1097-1770)

Reproducibility materials for the preprint:

> *Fibonacci Character Sums at Rank Index Two and Jacobsthal Sums*

> **Proof versus computation.** Mathematical results are proved in the manuscript.
> Repository computations provide finite-range verification of stated identities;
> they are not a substitute for the proofs.

---

## Citation and DOI records

| Object | Role | DOI |
| --- | --- | --- |
| Manuscript preprint | Theorems, proofs | [10.5281/zenodo.23264456](https://doi.org/10.5281/zenodo.23264456) |
| Software archive | Verifier + machine-readable evidence | [10.5281/zenodo.23263041](https://doi.org/10.5281/zenodo.23263041) |

Software record is maintained **manually** on Zenodo (New version per release).
A new Zenodo software version will be created manually from Git tag `v1.0.2` after the GitHub release source is finalized.

- Software metadata: [`CITATION.cff`](CITATION.cff), [`.zenodo.json`](.zenodo.json)
- Manuscript metadata: [`Paper/CITATION.cff`](Paper/CITATION.cff)
- Release notes: [`CHANGELOG.md`](CHANGELOG.md)

---

## Quick start

```bash
git clone https://github.com/Majid-Ghandali/fibonacci-character-sums-index-two.git
cd fibonacci-character-sums-index-two
python -m pip install -r Requirements.txt
```

Smoke test (does not replace the official full run; `--safe` writes to
`results/_smoke` and can never touch the official evidence directory):

```bash
python code/ghandali_verify_index2.py --safe --fresh
```

Equivalent explicit form:

```bash
python code/ghandali_verify_index2.py \
  --limit 10000 --sum-limit 10000 --jac-limit 10000 \
  --workers 1 --sum-indices 2,4 \
  --out results/reproduced-smoke-N10000 --fresh
```

Success ends with:

```text
script: ghandali_verify_index2.py 1.0.2
VERDICT: PASS
theorem mismatches: 0
```

Windows PowerShell (Python 3.12.x with numpy installed):

```powershell
python code\ghandali_verify_index2.py --safe --fresh
```

---

## Official evidence (version policy)

| Version | Evidence directory | Status |
| --- | --- | --- |
| 1.0.1 | `results/release-v1.0.1/` and `results/ghandali_index2/` | Historical, do not overwrite. Provenance: `release-v1.0.1` was produced by script 1.0.0 (see its log/summary); `ghandali_index2` is a later 1.0.1 rerun. See `CHANGELOG.md`. |
| **1.0.2** | `results/ghandali_index2_v1.0.2/` | **Official evidence.** Full run $N=2\cdot10^6$ completed 2026-10-13 (5036 s, Python 3.12.7, numpy 1.26.4); acceptance gate 53/53 passed. |

Do **not** rewrite 1.0.1 evidence. Mathematics of Theorems A and B is unchanged; 1.0.2 regenerates evidence under clean provenance (log version, theorem labels, `sha256sums.txt` written after the final log close).

### Full verification command (official 1.0.2)

Start from an empty `results\ghandali_index2_v1.0.2\` (the verifier refuses to
reuse a directory that already holds a finished run, and refuses N < 2,000,000):

```powershell
python code\ghandali_verify_index2.py --limit 2000000 --sum-limit 2000000 --jac-limit 2000000 --workers 14 --sum-indices 2,4 --out results\ghandali_index2_v1.0.2 --fresh
```

Gate after completion:

1. Log line: `script: ghandali_verify_index2.py 1.0.2`
2. `VERDICT: PASS`, `theorem_mismatch_total = 0`
3. Byte-match of `tested` counts against the 1.0.1 full-run summary for shared theorem keys (mathematics unchanged)
4. `sha256sums.txt` validates (it also covers the copy of the verifier stored in the directory):
   `cd results/ghandali_index2_v1.0.2 && sha256sum -c sha256sums.txt`

All four gates (and more) are automated:

```powershell
python code\gate_official.py --dir results\ghandali_index2_v1.0.2 --ref results\ghandali_index2\verification-summary.json --script code\ghandali_verify_index2.py
```

Headline counts at $N=2\cdot10^6$, confirmed by the official 1.0.2 run (2026-10-13; Python 3.12.7, numpy 1.26.4; runtime 5036 s; acceptance gate 53/53):

| Quantity | Value |
| --- | ---: |
| Rank inventory primes | 148,931 |
| Direct $S(p)$ records | 52,477 |
| Split $i_p=2$, $p\equiv1\pmod8$ | 6,551 |
| Split $i_p=2$, $p\equiv5\pmod8$ | 6,583 |
| Inert $i_p=2$, $p\equiv1\pmod8$ | 14,591 |
| Inert $i_p=2$, $p\equiv5\pmod8$ | 14,696 |
| Split $i_p=4$, $p\equiv5\pmod8$ | 6,697 |
| Split $i_p=4$, $p\equiv1\pmod8$ | 3,359 |
| Lehmer 1958 check (split, $i_p\in\{2,4\}$) | 23,190 |

SHA-256 of the verifier (also stored inside the evidence directory):
`7a74db81653b6f4fe5167d7f9d8639e9da2ad92d6cc58d1368f4047abaa9ef96`

---

## Verification scope

| Layer | Role |
| --- | --- |
| Theorem-level | Finite-range checks of identities proved in the manuscript (including signed inert law, index-four vanishing/parity) |
| Exploratory | Observations not claimed as theorems |
| Density context | Numerical ratios only; inert $10A/19$ is GRH-conditional (Roskam--Chen), not proved by the finite run |

---

## Repository layout

```text
Paper/                         manuscript (CC BY 4.0)
code/ghandali_verify_index2.py verifier (MIT), version 1.0.2
code/gate_official.py          post-run acceptance gate
results/ghandali_index2_v1.0.2/  official evidence after full run
CITATION.cff                   software citation
.zenodo.json                   Zenodo metadata (manual deposit)
CHANGELOG.md                   version history
Requirements.txt               numpy (the only third-party dependency)
LICENSE                        MIT (code/results)
Paper/LICENSE                  CC BY 4.0 (manuscript)
```

---

## Manuscript compilation

```bash
cd Paper
pdflatex -interaction=nonstopmode ghandali_index2.tex
bibtex ghandali_index2
pdflatex -interaction=nonstopmode ghandali_index2.tex
pdflatex -interaction=nonstopmode ghandali_index2.tex
```

Or: `.\compile.ps1` from the repository root.

---

## License

| Content | License |
| --- | --- |
| Code, scripts, verification outputs | MIT |
| Manuscript text, bibliography, PDF | CC BY 4.0 |

See [`LICENSES.md`](LICENSES.md).

---

## Author

**Majid Ghandali**
Independent Researcher, Tehran, Iran
Email: majid.ghandali@gmail.com
ORCID: [0009-0001-1097-1770](https://orcid.org/0009-0001-1097-1770)
