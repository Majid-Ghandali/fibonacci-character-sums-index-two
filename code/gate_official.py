"""
gate_official.py -- post-run acceptance gate for an official verifier run.

Usage (from the repository root):
  python code/gate_official.py --dir results/ghandali_index2_v1.0.2 \
      --ref results/ghandali_index2/verification-summary.json \
      --script code/ghandali_verify_index2.py

Exit code 0 = every gate passed, 1 = at least one gate failed.
Standard library only.
"""
import argparse, hashlib, json, sys
from pathlib import Path

REQUIRED = ["ghandali_index2_records.csv", "ghandali_verify_index2.py",
            "numerical-table.csv", "numerical-table.tex",
            "rank-index-all-primes.csv", "verification-run.log",
            "verification-summary.json", "verification-workbook.xlsx",
            "sha256sums.txt"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--ref", default=None, help="1.0.1 full-run verification-summary.json")
    ap.add_argument("--script", default=None, help="repo copy of the verifier to compare against")
    ap.add_argument("--expect-version", default="1.0.2")
    ap.add_argument("--expect-limit", type=int, default=2_000_000)
    a = ap.parse_args()
    d = Path(a.dir)
    fails, oks = [], []
    def ok(m): oks.append(m); print("  PASS  " + m)
    def bad(m): fails.append(m); print("  FAIL  " + m)
    def chk(c, m): ok(m) if c else bad(m)

    print(f"gate for {d}")
    chk(d.is_dir(), "output directory exists")
    if not d.is_dir():
        sys.exit(1)
    for f in REQUIRED:
        chk((d / f).is_file(), f"file present: {f}")
    chk(not (d / "archive").exists(), "no archive/ folder (directory was clean before the run)")
    chk(not list(d.glob("*.pkl")) and not list(d.glob("*.pkl.tmp")), "no leftover checkpoint")

    # manifest
    man = d / "sha256sums.txt"
    if man.is_file():
        lines = [l for l in man.read_text(encoding="utf-8").splitlines() if l.strip()]
        chk(len(lines) >= 8, f"manifest non-empty ({len(lines)} lines)")
        listed = set()
        for l in lines:
            h, _, name = l.partition("  ")
            listed.add(name)
            p = d / name
            chk(p.is_file() and sha(p) == h.lower(), f"manifest hash matches: {name}")
        actual = {p.relative_to(d).as_posix() for p in d.rglob("*") if p.is_file()} - {"sha256sums.txt"}
        chk(actual == listed, "manifest lists exactly the files present")

    # summary
    sp = d / "verification-summary.json"
    s = json.loads(sp.read_text(encoding="utf-8")) if sp.is_file() else {}
    m = s.get("metadata", {})
    chk(m.get("version") == a.expect_version, f"summary version == {a.expect_version} (got {m.get('version')})")
    chk(s.get("verdict") == "PASS", "verdict PASS")
    chk(s.get("theorem_mismatch_total") == 0, "theorem mismatches == 0")
    chk(s.get("exploratory_mismatch_total") == 0, "exploratory mismatches == 0")
    chk(m.get("limit") == m.get("sum_limit") == m.get("jac_limit") == a.expect_limit,
        f"limit == sum_limit == jac_limit == {a.expect_limit}")
    chk(m.get("sums_enabled") is True and m.get("jacobsthal_enabled") is True, "sums and Jacobsthal enabled")
    chk(m.get("sum_indices") == [2, 4], "sum_indices == [2, 4]")

    # log
    lp = d / "verification-run.log"
    log = lp.read_text(encoding="utf-8") if lp.is_file() else ""
    chk(f"script: ghandali_verify_index2.py {a.expect_version}" in log, "log names the expected script version")
    chk("--safe" not in log.split("started:")[0], "argv does not contain --safe")
    chk("VERDICT: PASS" in log, "log ends with VERDICT: PASS")
    sc = d / "ghandali_verify_index2.py"
    if sc.is_file():
        chk(f"script_sha256: {sha(sc)}" in log, "script_sha256 in log == hash of script copy in the directory")
        if a.script:
            chk(sha(Path(a.script)) == sha(sc), "script copy == repository script")

    # reference comparison (shared theorem keys only)
    if a.ref:
        r = json.loads(Path(a.ref).read_text(encoding="utf-8"))
        rt, nt = r["theorem_checks"], s.get("theorem_checks", {})
        shared = [k for k in rt if k in nt]
        chk(len(shared) >= 14, f"{len(shared)} shared theorem keys with the reference run")
        for k in shared:
            same = (rt[k]["tested"], rt[k]["matches"], rt[k]["mismatch"]) == (nt[k]["tested"], nt[k]["matches"], nt[k]["mismatch"])
            chk(same, f"tested/matches/mismatch identical to reference: {k} ({nt[k]['tested']})")
        rc = r["density_context_not_theorem"]["index2_counts_by_p_mod_8"]
        nc = s.get("density_context_not_theorem", {}).get("index2_counts_by_p_mod_8")
        chk(rc == nc, "index-two counts by p mod 8 identical to reference")
        chk(r["density_context_not_theorem"]["primes_counted"] == s.get("density_context_not_theorem", {}).get("primes_counted"),
            "rank-inventory prime count identical to reference")
        print("  info  new keys (not in reference): " + ", ".join(k for k in nt if k not in rt))

    print(f"\n{len(oks)} passed, {len(fails)} failed")
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
