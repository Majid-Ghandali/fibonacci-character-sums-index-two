# Changelog

## 1.0.2 (2026-10-13)
- Verifier hardened; **no mathematical check changed** (CSV/TeX data files are byte-identical to
  those of the pre-hardening 1.0.2 build for the same arguments).
- `--safe` without `--out` writes to `results/_smoke`; `--safe` can never target
  `results/ghandali_index2_v1.0.2/`.
- The official directory refuses N < 2,000,000 and refuses reuse once it holds a
  finished run (`--force-official` overrides).
- The verifier copies itself into the output directory; `sha256sums.txt` therefore
  covers the exact script that produced the evidence (`cd <dir>; sha256sum -c sha256sums.txt`).
- `sha256sums.txt` is no longer written empty.
- New theorem-level checks: signed inert law (S == -a* mod 5), index-four vanishing
  (Cor. 3.5) and parity (Lemma 3.1), E. Lehmer 1958 criterion.
- Evidence directory renamed: `results/ghandali_index2_v1.0.2/`; CSV name is now
  `ghandali_index2_records.csv`.
- Requirements reduced to numpy.
- Official full run (N = 2,000,000; `--workers 14`; 5036 s; Python 3.12.7, numpy 1.26.4):
  `VERDICT: PASS`, 0 theorem and 0 exploratory mismatches; acceptance gate
  `code/gate_official.py` passed 53/53 against the 1.0.1 rerun summary
  (verifier SHA-256 `7a74db81653b6f4fe5167d7f9d8639e9da2ad92d6cc58d1368f4047abaa9ef96`).
- Added `code/gate_official.py` (post-run acceptance gate).
- Root `LICENSE` is now the MIT text (it previously contained CC BY 4.0 text and
  referred to a non-existent `scripts/` directory); the manuscript stays CC BY 4.0
  in `Paper/LICENSE`.
- `.gitignore` now tracks `results/ghandali_index2_v1.0.2/` and ignores run
  leftovers (`archive/`, checkpoints).
- Manuscript: corrected the primitive-root argument in Proposition 4.1, relabelled
  the index-four parity statement as a lemma-level check, and removed the
  unsupported index-three remark.

## 1.0.1 (2026-10-09 / 2026-10-10)
- Packaging release. Provenance note: the files in `results/release-v1.0.1/` were
  produced by verifier script version 1.0.0 (see `verification-run.log` and
  `verification-summary.json` there), and its table lists the signed inert law as
  exploratory. `results/ghandali_index2/` is a later rerun with version 1.0.1.
  Both are kept unchanged as historical evidence; the official evidence for the
  manuscript is the 1.0.2 bundle.
- The Git tag `v1.0.1` contains a damaged `Paper/ghandali_index2.bib` (two entries
  without headers); fixed after the tag.
