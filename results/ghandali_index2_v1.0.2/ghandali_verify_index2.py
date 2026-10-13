# ghandali_verify_index2.py
"""
ghandali_verify_index2.py  (1.0.2)

Reproducibility verifier for:

    Fibonacci Character Sums at Rank Index Two and Jacobsthal Sums
    Majid Ghandali

Main corrections and upgrades
-----------------------------
1. Separates THEOREM checks from EXPLORATORY checks.
   (Signed inert law S == -a* (mod 5) is a THEOREM check.)
2. Index-four vanishing (split, p==5 mod 8) and index-four parity
   (split, p==1 mod 8) are THEOREM checks (Cor. 3.5 / Lemma 3.1).
3. Density context keys distinguish conjectural split weights from the
   GRH-conditional inert Dirichlet density 10A/19.
4. Logs python/numpy/platform, script sha-256, and argv; writes sha256sums.txt.
5. Checkpoints bind to script_version + script_sha256.
6. Progress bar, interactive UI, resume, atomic checkpoints.
7. 1.0.2: single-line progress bar clipped to terminal width (no
   wrapping); pipe-safe progress (one line per 10 %); checkpoint_format
   f-string fix; summary key 'jacobsthal_enabled'; interactive default
   workers = logical CPUs - 2; write_manifest no longer leaves sha256sums
   empty (append was unreachable under continue); SAFE mode writes to
   results/_smoke so official results/ghandali_index2_v1.0.2 is reserved
   for the controlled full run. Mathematics and checks unchanged.
8. 1.0.2 (final hardening, no change to any mathematical check):
   --safe without --out now writes to results/_smoke (previously it wrote
   to the official directory); --safe can never target the official
   directory; the official directory refuses to be reused once it holds a
   finished run (verification-summary.json) and refuses runs with
   N < 2,000,000 unless --force-official is given; the verifier copies
   itself into the output directory so sha256sums.txt also covers the exact
   script that produced the evidence; interactive modes no longer default
   to the official directory for N < 2,000,000.

Usage
-----
Interactive UI:

    python ghandali_verify_index2.py

Safe smoke test:

    python ghandali_verify_index2.py --safe

Moderate verification:

    python ghandali_verify_index2.py --limit 1000000 --sum-limit 1000000 --jac-limit 200000

Fresh run, ignoring checkpoint:

    python ghandali_verify_index2.py --limit 1000000 --fresh

Resume automatically if checkpoint exists:

    python ghandali_verify_index2.py --limit 1000000 --resume

Requirements
------------
Python 3.12.7 (tested) with numpy>=1.26.4,<2.

License
-------
MIT.
"""

import os

for _v in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_v, "1")

import zipfile
import time
import argparse
import csv
import datetime as _dt
import hashlib
import json
import math
import platform
import pickle
import shutil
import sys
from multiprocessing import Pool, cpu_count, freeze_support
from pathlib import Path

import numpy as np


VERSION = "1.0.2"
OFFICIAL_OUT = f"results/ghandali_index2_v{VERSION}"
SMOKE_OUT = "results/_smoke"
OFFICIAL_MIN_LIMIT = 2_000_000
TITLE = "Fibonacci Character Sums at Rank Index Two and Jacobsthal Sums"
AUTHOR = "Majid Ghandali"

A_ARTIN = 0.3739558136192022880547280543464164151116292486061500


# ============================================================
# Basic arithmetic
# ============================================================

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    v = pow(a, (p - 1) // 2, p)
    if v == 1:
        return 1
    if v == p - 1:
        return -1
    raise ArithmeticError(f"unexpected Legendre pow value for a={a}, p={p}: {v}")


def quartic_symbol(a, p):
    if p % 4 != 1:
        return None
    a %= p
    if a == 0:
        return 0
    if legendre(a, p) != 1:
        return None
    v = pow(a, (p - 1) // 4, p)
    if v == 1:
        return 1
    if v == p - 1:
        return -1
    return None


def delta_split(p):
    r = p % 5
    if r == 1:
        return 1
    if r == 4:
        return -1
    return None


def fib_pair(n, p):
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1, p)
    c = (a * ((2 * b - a) % p)) % p
    d = (a * a + b * b) % p
    if n & 1:
        return d, (c + d) % p
    return c, d


def fib_mod(n, p):
    return fib_pair(n, p)[0]


def primes_upto(n):
    if n < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    r = math.isqrt(n)
    for q in range(2, r + 1):
        if sieve[q]:
            sieve[q * q:n + 1:q] = False
    return np.nonzero(sieve)[0].astype(np.int64)


def factor_distinct(n, small_primes):
    m = int(n)
    out = []
    for q in small_primes:
        q = int(q)
        if q * q > m:
            break
        if m % q == 0:
            out.append(q)
            while m % q == 0:
                m //= q
    if m > 1:
        out.append(m)
    return out


def rank_of_apparition(p, small_primes):
    eps5 = legendre(5, p)
    n = p - eps5
    d = n
    for q in factor_distinct(n, small_primes):
        while d % q == 0 and fib_mod(d // q, p) == 0:
            d //= q
    return d, n // d, eps5


def rank_minimality_check(p, z, n, small_primes):
    if n % z != 0:
        return False
    if fib_mod(z, p) != 0:
        return False
    for q in factor_distinct(z, small_primes):
        if z % q == 0 and fib_mod(z // q, p) == 0:
            return False
    return True


def two_square_data(p):
    if p % 4 != 1:
        return None, None
    r = math.isqrt(p)
    for a_abs in range(1, r + 1, 2):
        b2 = p - a_abs * a_abs
        if b2 < 0:
            break
        b = math.isqrt(b2)
        if b > 0 and b * b == b2 and b % 2 == 0:
            a_star = a_abs if a_abs % 4 == 1 else -a_abs
            return a_star, b
    raise ArithmeticError(f"two-square representation not found for p={p}")


# ============================================================
# Direct Sums
# ============================================================

def chi_table_numpy(p):
    chi = np.full(p, -1, dtype=np.int8)
    chi[0] = 0
    xs = np.arange(1, (p + 1) // 2, dtype=np.int64)
    chi[(xs * xs) % p] = 1
    return chi


def fibonacci_character_sum(p, z, block=131072):
    chi = chi_table_numpy(p)

    B = int(min(block, z))
    base_f = np.empty(B, dtype=np.int64)
    base_fp1 = np.empty(B, dtype=np.int64)

    fj, fj1 = 1 % p, 1 % p
    for k in range(B):
        base_f[k] = fj
        base_fp1[k] = fj1
        fj, fj1 = fj1, (fj + fj1) % p

    total = 0
    m = 0
    Fm = 0
    Fm_minus1 = 1

    while m < z:
        L = min(B, z - m)
        vals = (Fm * base_fp1[:L] + Fm_minus1 * base_f[:L]) % p
        total += int(chi[vals].sum(dtype=np.int64))

        m += L
        Fm, Fm_plus1 = fib_pair(m, p)
        Fm_minus1 = (Fm_plus1 - Fm) % p

    return int(total)


def jacobsthal_sum(p, c, chunk=1_000_000):
    chi = chi_table_numpy(p)
    c %= p
    total = 0
    for start in range(0, p, chunk):
        stop = min(p, start + chunk)
        x = np.arange(start, stop, dtype=np.int64)
        y = (x * x + c) % p
        prod = chi[x].astype(np.int16) * chi[y].astype(np.int16)
        total += int(prod.sum(dtype=np.int64))
    return int(total)



# ============================================================
# Progress bar
# ============================================================

_ACTIVE_BAR = None


class ProgressBar:

    def __init__(self, label, total, width=42, enabled=True):
        global _ACTIVE_BAR
        self.label = str(label)
        self.total = int(total)
        self.max_width = int(width)
        self.enabled = bool(enabled and self.total >= 0)
        self.start = time.time()
        self.last_draw = 0.0
        self.finished = False
        self.on_screen = False
        self.last_decile = -1
        try:
            self.is_tty = bool(sys.stdout.isatty())
        except Exception:
            self.is_tty = False

        enc = (getattr(sys.stdout, "encoding", None) or "").lower()
        self.utf = "utf" in enc
        if self.utf:
            self.fill, self.empty = "\u2588", "\u2591"
        else:
            self.fill, self.empty = "#", "-"

        if self.enabled:
            _ACTIVE_BAR = self

    @staticmethod
    def _columns():
        try:
            return max(20, shutil.get_terminal_size((80, 24)).columns)
        except Exception:
            return 80

    def _compose(self, done, note, now, cols):
        frac = 1.0 if self.total <= 0 else max(0.0, min(1.0, done / self.total))
        pct = 100.0 * frac
        elapsed = now - self.start
        rate = f"{done / elapsed:.1f}/s" if done > 0 and elapsed > 0 else "--/s"

        eta = ""
        if 0 < done < self.total and elapsed > 0:
            m, s = divmod(int((self.total - done) * elapsed / done), 60)
            h, m = divmod(m, 60)
            eta = f"  ETA {h:d}:{m:02d}:{s:02d}"

        note = str(note or "")
        if not self.utf:
            note = note.replace("\u2264", "<=").encode("ascii", "replace").decode()

        tw = len(str(self.total))
        left = f"{self.label:<14} |"
        right = f"| {pct:6.2f}%  {done:>{tw}d}/{self.total}  {rate}{eta}"
        if note:
            right += "  " + note

        avail = cols - 1                     # never write into the last column
        bar_w = min(self.max_width, avail - len(left) - len(right))
        if bar_w < 10:                       # narrow console: drop extras
            right = f"| {pct:6.2f}%  {done}/{self.total}"
            bar_w = min(self.max_width, avail - len(left) - len(right))
        bar_w = max(0, bar_w)

        filled = int(round(bar_w * frac))
        bar = self.fill * filled + self.empty * (bar_w - filled)
        return (left + bar + right)[:avail], pct

    def update(self, done, note="", force=False):
        if not self.enabled or self.finished:
            return
        done = int(done)
        now = time.time()
        if not force and now - self.last_draw < 0.15 and done < self.total:
            return
        self.last_draw = now

        cols = self._columns()
        line, pct = self._compose(done, note, now, cols)

        if self.is_tty:
            sys.stdout.write("\r" + line.ljust(cols - 1))
            self.on_screen = True
        else:
            decile = int(pct // 10)
            if decile == self.last_decile and not force:
                return
            self.last_decile = decile
            sys.stdout.write(line + "\n")
        sys.stdout.flush()

    def clear(self):
        """Erase the bar so that an ordinary print() does not mix with it."""
        if self.enabled and self.is_tty and self.on_screen and not self.finished:
            sys.stdout.write("\r" + " " * (self._columns() - 1) + "\r")
            sys.stdout.flush()
            self.on_screen = False

    def finish(self, note="done"):
        global _ACTIVE_BAR
        if not self.enabled or self.finished:
            return
        self.update(self.total, note=note, force=True)
        if self.is_tty:
            sys.stdout.write("\n")
            sys.stdout.flush()
        self.finished = True
        if _ACTIVE_BAR is self:
            _ACTIVE_BAR = None


# ============================================================
# Logging and files
# ============================================================

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def write_manifest(out_dir):
    out_dir = Path(out_dir)
    manifest = out_dir / "sha256sums.txt"
    lines = []
    for p in sorted(out_dir.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(out_dir).as_posix()
        if (
            rel == "sha256sums.txt"
            or rel.startswith("archive/")
            or rel.endswith(".pkl")
            or rel.endswith(".pkl.tmp")
        ):
            continue
        lines.append("%s  %s" % (sha256_file(p).lower(), rel))
    with open(manifest, "w", encoding="utf-8", newline=chr(10)) as fh:
        fh.write(chr(10).join(lines) + chr(10))
    return manifest



class Logger:
    def __init__(self, path):
        self.path = Path(path)
        self.fh = open(self.path, "w", encoding="utf-8", newline="\n")

    def __call__(self, msg):
        line = str(msg)
        if _ACTIVE_BAR is not None:
            _ACTIVE_BAR.clear()
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()

    def close(self):
        self.fh.close()


def archive_existing_outputs(out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    suffixes = {".json", ".csv", ".tex", ".log", ".xlsx", ".txt"}
    existing = [
        p for p in out_dir.iterdir()
        if p.is_file() and p.suffix.lower() in suffixes
    ]
    if not existing:
        return None

    stamp = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    arch = out_dir / "archive" / stamp
    k = 0
    while arch.exists():
        k += 1
        arch = out_dir / "archive" / f"{stamp}-{k}"
    arch.mkdir(parents=True, exist_ok=False)

    for p in existing:
        shutil.move(str(p), str(arch / p.name))
    return arch


# ============================================================
# Checkpoint
# ============================================================

def checkpoint_path(out_dir):
    # Checkpoint is a trusted local file only (not for network exchange).
    return Path(out_dir) / "checkpoint_ghandali_index2.pkl"


def make_checkpoint_meta(args, sum_indices):
    try:
        _script_sha = hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest()
    except OSError:
        _script_sha = "unavailable"
    return {
        "checkpoint_format": f"ghandali_index2_v{VERSION}",
        "script_version": VERSION,
        "script_sha256": _script_sha,
        "limit": int(args.limit),
        "sum_limit": int(args.sum_limit),
        "jac_limit": int(args.jac_limit),
        "sum_indices": sorted(int(x) for x in sum_indices),
        "no_sums": bool(args.no_sums),
        "no_jac": bool(args.no_jac),
        "rank_chunk": int(args.rank_chunk),
        "block": int(args.block),
        "jac_chunk": int(args.jac_chunk),
    }


def save_checkpoint(out_dir, meta, rank_parts, recs, phase):
    path = checkpoint_path(out_dir)
    payload = {
        "meta": meta,
        "phase": str(phase),
        "rank_parts": rank_parts,
        "recs": recs,
        "saved_at": _dt.datetime.now().isoformat(timespec="seconds"),
    }

    tmp = path.with_suffix(".pkl.tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
        fh.flush()
        os.fsync(fh.fileno())

    try:
        os.replace(str(tmp), str(path))
    except OSError:
        if path.exists():
            path.unlink()
        shutil.move(str(tmp), str(path))


def load_checkpoint(out_dir, meta):
    path = checkpoint_path(out_dir)
    if not path.exists():
        return None, "missing"

    try:
        with open(path, "rb") as fh:
            payload = pickle.load(fh)
    except Exception as exc:
        return None, f"unreadable: {exc}"

    if payload.get("meta") != meta:
        return None, "meta-mismatch"

    return payload, "ok"


def clear_checkpoint(out_dir):
    path = checkpoint_path(out_dir)
    if path.exists():
        path.unlink()


def ask_resume_or_fresh(payload):
    phase = payload.get("phase", "unknown")
    rank_parts = payload.get("rank_parts", {})
    recs = payload.get("recs", [])
    saved_at = payload.get("saved_at", "unknown")

    print()
    print("Existing checkpoint found.")
    print(f"  saved_at:   {saved_at}")
    print(f"  phase:      {phase}")
    print(f"  rank parts: {len(rank_parts)}")
    print(f"  S records:  {len(recs)}")
    print()
    ans = input("Resume previous run? [Y/n]: ").strip().lower()
    return ans in ("", "y", "yes")


# ============================================================
# Interactive UI
# ============================================================

def _cpu_default():
    return max(1, (os.cpu_count() or 3) - 2)


def _ask_int(prompt, default, minimum=None, maximum=None):
    while True:
        raw = input(f"{prompt} [{default}]: ").strip()
        if raw == "":
            val = int(default)
        else:
            try:
                val = int(raw.replace("_", ""))
            except ValueError:
                print("  invalid integer; try again.")
                continue

        if minimum is not None and val < minimum:
            print(f"  must be >= {minimum}")
            continue
        if maximum is not None and val > maximum:
            print(f"  must be <= {maximum}")
            continue
        return val


def _ask_choice(prompt, choices, default):
    choices = list(choices)
    while True:
        raw = input(f"{prompt} [{default}]: ").strip().lower()
        if raw == "":
            raw = default
        if raw in choices:
            return raw
        print("  choices: " + ", ".join(choices))


def interactive_setup_argv():
    print()
    print("=" * 72)
    print(" Fibonacci Character Sums at Rank Index Two -- verifier UI")
    print("=" * 72)
    print()
    print("Modes:")
    print("  1) SAFE smoke test")
    print("  2) Standard verification")
    print("  3) Rank/density inventory only -- no direct S Sums")
    print("  4) Custom")
    print()

    mode = _ask_choice("Select mode", ["1", "2", "3", "4"], "2")
    argv = []

    out_default = OFFICIAL_OUT

    if mode == "1":
        # Smoke must never land in the official versioned evidence directory.
        argv += ["--safe", "--out", "results/_smoke"]
        return argv

    if mode == "2":
        limit = _ask_int("Prime limit N", 1_000_000, minimum=100)
        sum_limit = _ask_int("Direct S(p) limit", limit, minimum=0, maximum=limit)
        jac_limit = _ask_int("Jacobsthal direct-check limit", min(200_000, limit), minimum=0, maximum=limit)

        print()
        print("Worker guidance:")
        print("  workers=1 is safest on laptops; default = auto (logical CPUs - 2).")
        print(f"  auto suggestion for this system: {_cpu_default()}")
        workers = _ask_int(
            "Worker processes for rank inventory",
            _cpu_default(),
            minimum=1,
            maximum=max(1, (os.cpu_count() or 4) * 2),
        )

        indices_raw = input("Rank indices to compute direct S for [2,4]: ").strip()
        indices = indices_raw if indices_raw else "2,4"

        out_def = out_default if limit >= OFFICIAL_MIN_LIMIT else f"results/reproduced-N{limit}"
        out = input(f"Output directory [{out_def}]: ").strip() or out_def

        argv += [
            "--limit", str(limit),
            "--sum-limit", str(sum_limit),
            "--jac-limit", str(jac_limit),
            "--workers", str(workers),
            "--sum-indices", indices,
            "--out", out,
        ]
        return argv

    if mode == "3":
        limit = _ask_int("Prime limit N", 10_000_000, minimum=100)
        workers = _ask_int(
            "Worker processes for rank inventory",
            _cpu_default(),
            minimum=1,
            maximum=max(1, (os.cpu_count() or 4) * 2),
        )
        out_def = out_default if limit >= OFFICIAL_MIN_LIMIT else f"results/reproduced-N{limit}"
        out = input(f"Output directory [{out_def}]: ").strip() or out_def

        argv += [
            "--limit", str(limit),
            "--workers", str(workers),
            "--no-sums",
            "--out", out,
        ]
        return argv

    # mode 4 custom
    limit = _ask_int("Prime limit N", 1_000_000, minimum=100)
    sum_limit = _ask_int("Direct S(p) limit", limit, minimum=0, maximum=limit)
    jac_limit = _ask_int("Jacobsthal direct-check limit", min(200_000, limit), minimum=0, maximum=limit)
    workers = _ask_int("Worker processes", _cpu_default(), minimum=1, maximum=max(1, (os.cpu_count() or 3) - 2))
    rank_chunk = _ask_int("Rank checkpoint chunk size", 5000, minimum=100)
    ckpt_every = _ask_int("Checkpoint every K units", 25, minimum=1)
    indices_raw = input("Rank indices to compute direct S for [2,4]: ").strip()
    indices = indices_raw if indices_raw else "2,4"
    out_def = out_default if limit >= OFFICIAL_MIN_LIMIT else f"results/reproduced-N{limit}"
    out = input(f"Output directory [{out_def}]: ").strip() or out_def

    no_jac = _ask_choice("Disable Jacobsthal Sums?", ["y", "n"], "n")
    no_progress = _ask_choice("Disable progress bars?", ["y", "n"], "n")

    argv += [
        "--limit", str(limit),
        "--sum-limit", str(sum_limit),
        "--jac-limit", str(jac_limit),
        "--workers", str(workers),
        "--rank-chunk", str(rank_chunk),
        "--checkpoint-every", str(ckpt_every),
        "--sum-indices", indices,
        "--out", out,
    ]
    if no_jac == "y":
        argv.append("--no-jac")
    if no_progress == "y":
        argv.append("--no-progress")
    return argv


# ============================================================
# Worker
# ============================================================

def _rank_worker(payload):
    idx, prime_chunk, small_primes = payload
    rows = []
    for pp in prime_chunk:
        p = int(pp)
        if p < 3 or p == 5:
            continue
        z, ip, eps5 = rank_of_apparition(p, small_primes)

        if z % 2 == 1:
            chi_fzp1 = legendre(fib_mod(z + 1, p), p)
            chiok = 1 if chi_fzp1 == legendre(2, p) else 0
        else:
            chiok = -1

        rows.append((p, eps5, z, ip, chiok))
    return int(idx), rows


# ============================================================
# Prediction helpers
# ============================================================

def inert_sign_rule_prediction(a_star, b):
    target = (-a_star) % 5
    if b % 5 == target:
        return b
    if (-b) % 5 == target:
        return -b
    return None


def main_theorem_prediction_for_record(rec):
    p = rec["p"]
    eps5 = rec["eps5"]
    ip = rec["index"]
    S = rec.get("S")
    a_star = rec.get("a_star")
    b = rec.get("b")
    J5 = rec.get("J5")

    if S is None:
        return None, None

    if ip == 2 and eps5 == 1:
        d = delta_split(p)
        if p % 8 == 1:
            pred = d * a_star
            return pred, int(S == pred)
        if p % 8 == 5:
            pred = -d
            return pred, int(S == pred)

    if ip == 2 and eps5 == -1:
        if p % 8 == 1:
            spred = inert_sign_rule_prediction(a_star, b)
            if spred is None:
                return None, 0
            if J5 is not None:
                ok = (S == spred) and (2 * S == J5) and (abs(S) == b)
            else:
                ok = (S == spred) and (abs(S) == b)
            return spred, int(ok)
        if p % 8 == 5:
            pred = 0
            return pred, int(S == 0)

    return None, None


def build_sum_record(row, small_primes, args):
    p, eps5, z, ip, chiok = [int(x) for x in row]

    rec = {
        "p": p,
        "eps5": eps5,
        "z": z,
        "index": ip,
        "p_mod_8": p % 8,
        "rank_ok": int(rank_minimality_check(p, z, p - eps5, small_primes)),
        "chi_F_z_plus_1_check_odd_z": chiok,
        "S": None,
        "a_star": None,
        "b": None,
        "delta": delta_split(p) if eps5 == 1 else None,
        "main_theorem_prediction": None,
        "main_theorem_match": None,
        "J5": None,
        "Jm5": None,
        "J10": None,
        "Jm10": None,
        "Jm1": None,
        "inert_sign_rule_prediction": None,
        "inert_sign_rule_match": None,
        "classification": "",
    }

    if p % 4 == 1:
        a_star, b = two_square_data(p)
        rec["a_star"] = a_star
        rec["b"] = b

    S = fibonacci_character_sum(p, z, block=args.block)
    rec["S"] = S

    if ip == 2:
        rec["classification"] = "THEOREM-index-two"
    elif ip == 4 and eps5 == 1 and p % 8 == 5:
        rec["classification"] = "THEOREM-index-four-p5mod8"  # Cor. 3.5(b)
    elif ip == 4 and eps5 == 1 and p % 8 == 1:
        rec["classification"] = "THEOREM-parity-index-four-p1mod8"  # Lemma 3.1
    else:
        rec["classification"] = "computed"

    if (not args.no_jac) and p <= args.jac_limit:
        if ip == 2 and p % 4 == 1:
            rec["Jm1"] = jacobsthal_sum(p, -1, chunk=args.jac_chunk)

        if ip == 2 and eps5 == -1 and p % 8 == 1:
            rec["J5"] = jacobsthal_sum(p, 5, chunk=args.jac_chunk)
            rec["Jm5"] = jacobsthal_sum(p, -5, chunk=args.jac_chunk)
            rec["J10"] = jacobsthal_sum(p, 10, chunk=args.jac_chunk)
            rec["Jm10"] = jacobsthal_sum(p, -10, chunk=args.jac_chunk)

    pred, ok = main_theorem_prediction_for_record(rec)
    rec["main_theorem_prediction"] = pred
    rec["main_theorem_match"] = ok

    if ip == 2 and eps5 == -1 and p % 8 == 1 and rec["a_star"] is not None:
        spred = inert_sign_rule_prediction(rec["a_star"], rec["b"])
        rec["inert_sign_rule_prediction"] = spred
        rec["inert_sign_rule_match"] = None if spred is None else int(S == spred)

    return rec


# ============================================================
# Aggregation
# ============================================================

def make_check(tested, bad):
    bad = list(bad)
    return {
        "tested": int(tested),
        "matches": int(tested - len(bad)),
        "mismatch": int(len(bad)),
        "first_mismatches": [int(x) for x in bad[:10]],
    }


def aggregate_checks(R, recs):
    theorem = {}
    exploratory = {}

    if len(R) == 0:
        return {}, theorem, exploratory, 0.0

    p = R[:, 0].astype(np.int64)
    eps5 = R[:, 1].astype(np.int64)
    z = R[:, 2].astype(np.int64)
    ip = R[:, 3].astype(np.int64)
    chiok = R[:, 4].astype(np.int64)

    P = int(len(p))

    bad_integrality = p[(z * ip) != (p - eps5)]
    theorem["rank_index_integrality_z_times_i_eq_n"] = make_check(P, bad_integrality)

    bad_parity = p[((ip % 2) == 0) != ((p % 4) == 1)]
    theorem["index_parity_i_even_iff_p_1mod4"] = make_check(P, bad_parity)

    odd_mask = chiok >= 0
    bad_chi = p[odd_mask & (chiok == 0)]
    theorem["odd_z_chi_F_z_plus_1_eq_2_over_p"] = make_check(int(odd_mask.sum()), bad_chi)

    qmask = (ip == 2) & (eps5 == 1) & ((p % 4) == 1)
    bad_q = []
    for pp in p[qmask]:
        pp = int(pp)
        lhs = quartic_symbol(5, pp)
        rhs = -delta_split(pp) * legendre(2, pp)
        if lhs != rhs:
            bad_q.append(pp)
    theorem["quartic_5_eq_minus_delta_times_2_over_p_split_idx2"] = make_check(
        int(qmask.sum()), bad_q
    )

    def pick(index=None, eps=None, m8=None, need_S=True, need=None, pmod4=None):
        out = []
        for r in recs:
            if index is not None and r["index"] != index:
                continue
            if eps is not None and r["eps5"] != eps:
                continue
            if m8 is not None and r["p"] % 8 != m8:
                continue
            if pmod4 is not None and r["p"] % 4 != pmod4:
                continue
            if need_S and r.get("S") is None:
                continue
            if need is not None and r.get(need) is None:
                continue
            out.append(r)
        return out

    def chk_rec(target, name, rows, pred):
        bad = []
        for r in rows:
            try:
                ok = bool(pred(r))
            except Exception:
                ok = False
            if not ok:
                bad.append(r["p"])
        target[name] = make_check(len(rows), bad)

    refl_rows = [r for r in recs if r.get("S") is not None and r["z"] % 2 == 1]
    chk_rec(
        theorem,
        "reflection_odd_z_p1mod4_and_S_eq_2overp_times_S",
        refl_rows,
        lambda r: (r["p"] % 4 == 1)
        and (r["S"] == legendre(2, r["p"]) * r["S"]),
    )

    refl_zero_rows = [
        r for r in recs
        if r.get("S") is not None and r["z"] % 2 == 1 and r["p"] % 8 == 5
    ]
    chk_rec(
        theorem,
        "reflection_odd_z_p5mod8_S_eq_0",
        refl_zero_rows,
        lambda r: r["S"] == 0,
    )

    chk_rec(
        theorem,
        "thmA_split_idx2_p1mod8_S_eq_delta_astar",
        pick(index=2, eps=1, m8=1),
        lambda r: r["S"] == delta_split(r["p"]) * r["a_star"],
    )

    chk_rec(
        theorem,
        "thmA_split_idx2_p5mod8_S_eq_minus_delta",
        pick(index=2, eps=1, m8=5),
        lambda r: r["S"] == -delta_split(r["p"]),
    )

    chk_rec(
        theorem,
        "jacobsthal_Jm1_eq_minus_2_times_2overp_times_astar_idx2_p1mod4",
        pick(index=2, pmod4=1, need="Jm1"),
        lambda r: r["Jm1"] == -2 * legendre(2, r["p"]) * r["a_star"],
    )

    chk_rec(
        theorem,
        "split_reduction_p1mod8_S_eq_quartic5_half_Jm1",
        pick(index=2, eps=1, m8=1, need="Jm1"),
        lambda r: (r["Jm1"] % 2 == 0)
        and (r["S"] == quartic_symbol(5, r["p"]) * (r["Jm1"] // 2)),
    )

    chk_rec(
        theorem,
        "split_reduction_p5mod8_S_eq_minus_quartic5",
        pick(index=2, eps=1, m8=5),
        lambda r: r["S"] == -quartic_symbol(5, r["p"]),
    )

    chk_rec(
        theorem,
        "thmB_inert_idx2_p1mod8_abs_S_eq_b",
        pick(index=2, eps=-1, m8=1),
        lambda r: abs(r["S"]) == r["b"],
    )

    chk_rec(
        theorem,
        "thmB_inert_idx2_p1mod8_2S_eq_J5",
        pick(index=2, eps=-1, m8=1, need="J5"),
        lambda r: 2 * r["S"] == r["J5"],
    )

    chk_rec(
        theorem,
        "thmB_inert_idx2_p1mod8_J5_eq_Jm5",
        pick(index=2, eps=-1, m8=1, need="Jm5"),
        lambda r: r["J5"] == r["Jm5"],
    )

    chk_rec(
        theorem,
        "thmB_inert_idx2_p5mod8_S_eq_0",
        pick(index=2, eps=-1, m8=5),
        lambda r: r["S"] == 0,
    )

    chk_rec(
        theorem,
        "thmB_inert_idx2_p1mod8_sign_rule_mod5",
        pick(index=2, eps=-1, m8=1),
        lambda r: r["inert_sign_rule_prediction"] is not None
        and r["S"] == r["inert_sign_rule_prediction"],
    )

    chk_rec(
        theorem,
        "lehmer1958_quartic5_iff_5_divides_b_split",
        pick(eps=1, pmod4=1),
        lambda r: r.get("b") is not None
        and (quartic_symbol(5, r["p"]) == 1) == (r["b"] % 5 == 0),
    )

    chk_rec(
        exploratory,
        "exploratory_inert_idx2_p1mod8_J10_eq_quartic2_J5",
        pick(index=2, eps=-1, m8=1, need="J10"),
        lambda r: r["J10"] == quartic_symbol(2, r["p"]) * r["J5"],
    )

    chk_rec(
        exploratory,
        "exploratory_inert_idx2_p1mod8_abs_Jpm10_eq_2absS",
        pick(index=2, eps=-1, m8=1, need="Jm10"),
        lambda r: abs(r["J10"]) == 2 * abs(r["S"]) == abs(r["Jm10"]),
    )

    chk_rec(
        theorem,
        "cor35_split_idx4_p5mod8_S_eq_0",
        pick(index=4, eps=1, m8=5),
        lambda r: r["S"] == 0,
    )

    idx4_p1_rows = pick(index=4, eps=1, m8=1)
    chk_rec(
        theorem,
        "lemma31_split_idx4_p1mod8_S_odd",
        idx4_p1_rows,
        lambda r: abs(r["S"]) % 2 == 1,
    )

    max_ratio = 0.0
    if idx4_p1_rows:
        max_ratio = max(abs(r["S"]) / math.sqrt(r["p"]) for r in idx4_p1_rows)

    exploratory["exploratory_split_idx4_p1mod8_max_absS_over_sqrtp"] = {
        "tested": len(idx4_p1_rows),
        "value": float(max_ratio),
        "claim_status": "exploratory_bound_observation_only",
    }

    def ratio(mask):
        return float(mask.sum()) / P / A_ARTIN if P else 0.0

    def by_mod(mask, mod):
        vals, cnt = np.unique(p[mask] % mod, return_counts=True)
        return {str(int(v)): int(c) for v, c in zip(vals, cnt)}

    vals, cnt = np.unique(ip, return_counts=True)
    count_by_index_first20 = {
        str(int(a)): int(b)
        for a, b in zip(vals, cnt)
        if int(a) <= 20
    }

    density = {
        "claim_status": "context_only_not_a_theorem_check",
        "primes_counted": P,
        "full_rank_over_A": ratio(ip == 1),
        "full_rank_split_over_A": ratio((ip == 1) & (eps5 == 1)),
        "full_rank_inert_over_A": ratio((ip == 1) & (eps5 == -1)),
        "context_pred_full_rank_split_over_A_CONDITIONAL": 9 / 19,
        "context_pred_full_rank_inert_over_A_CONDITIONAL": 10 / 19,
        "index2_split_over_A": ratio((ip == 2) & (eps5 == 1)),
        "index2_inert_over_A": ratio((ip == 2) & (eps5 == -1)),
        "context_pred_index2_split_over_A_CONJECTURAL": 9 / 38,
        "context_pred_index2_inert_over_A_DIRICHLET_CONDITIONAL": 10 / 19,
        "index2_counts_by_p_mod_8": {
            "split": by_mod((ip == 2) & (eps5 == 1), 8),
            "inert": by_mod((ip == 2) & (eps5 == -1), 8),
        },
        "index4_counts_by_p_mod_8": {
            "split": by_mod((ip == 4) & (eps5 == 1), 8),
            "inert": by_mod((ip == 4) & (eps5 == -1), 8),
        },
        "count_by_index_first20": count_by_index_first20,
    }

    return density, theorem, exploratory, max_ratio


# ============================================================
# Output writers
# ============================================================

def latex_escape(s):
    s = str(s)
    repl = {
        "_": r"\_",
        "%": r"\%",
        "&": r"\&",
        "#": r"\#",
        "{": r"\{",
        "}": r"\}",
    }
    for a, b in repl.items():
        s = s.replace(a, b)
    return s


def build_paper_table(theorem, exploratory):
    rows = []

    def add(category, name, statement, checks, status):
        c = checks.get(name, {"tested": 0, "matches": 0, "mismatch": 0})
        rows.append({
            "category": category,
            "name": name,
            "statement": statement,
            "tested": int(c.get("tested", 0)),
            "matches": int(c.get("matches", 0)),
            "mismatches": int(c.get("mismatch", 0)),
            "status": status,
        })

    add(
        "THEOREM",
        "thmA_split_idx2_p1mod8_S_eq_delta_astar",
        "split i_p=2, p=1 mod 8: S=delta_p a*",
        theorem,
        "proved-claim-check",
    )
    add(
        "THEOREM",
        "thmA_split_idx2_p5mod8_S_eq_minus_delta",
        "split i_p=2, p=5 mod 8: S=-delta_p",
        theorem,
        "proved-claim-check",
    )
    add(
        "THEOREM",
        "thmB_inert_idx2_p1mod8_abs_S_eq_b",
        "inert i_p=2, p=1 mod 8: |S|=b",
        theorem,
        "proved-claim-check",
    )
    add(
        "THEOREM",
        "thmB_inert_idx2_p1mod8_2S_eq_J5",
        "inert i_p=2, p=1 mod 8: 2S=Jac(5)",
        theorem,
        "proved-claim-check-where-J-computed",
    )
    add(
        "THEOREM",
        "thmB_inert_idx2_p5mod8_S_eq_0",
        "inert i_p=2, p=5 mod 8: S=0",
        theorem,
        "proved-claim-check",
    )
    add(
        "THEOREM",
        "thmB_inert_idx2_p1mod8_sign_rule_mod5",
        "inert i_p=2, p=1 mod 8: S in {+/-b}, S == -a* mod 5",
        theorem,
        "proved-claim-check",
    )
    add(
        "THEOREM",
        "cor35_split_idx4_p5mod8_S_eq_0",
        "split i_p=4, p=5 mod 8: S=0 (Cor. 3.5)",
        theorem,
        "proved-claim-check",
    )
    add(
        "THEOREM",
        "lemma31_split_idx4_p1mod8_S_odd",
        "split i_p=4, p=1 mod 8: S odd (Lemma 3.1)",
        theorem,
        "proved-claim-check",
    )
    add(
        "THEOREM",
        "lehmer1958_quartic5_iff_5_divides_b_split",
        "split: (5/p)_4 = 1 iff 5 | b (E. Lehmer 1958)",
        theorem,
        "proved-claim-check",
    )

    return rows


def write_rank_csv(R, path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("p,eps5,z,index,chi_F_z_plus_1_check_odd_z\n")
        np.savetxt(fh, R, fmt="%d", delimiter=",")


def write_records_csv(recs, path):
    fields = [
        "p",
        "eps5",
        "z",
        "index",
        "p_mod_8",
        "rank_ok",
        "chi_F_z_plus_1_check_odd_z",
        "S",
        "a_star",
        "b",
        "delta",
        "main_theorem_prediction",
        "main_theorem_match",
        "J5",
        "Jm5",
        "J10",
        "Jm10",
        "Jm1",
        "inert_sign_rule_prediction",
        "inert_sign_rule_match",
        "classification",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in recs:
            w.writerow(r)


def write_table_csv(rows, path):
    fields = ["category", "name", "statement", "tested", "matches", "mismatches", "status"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_table_tex(rows, path):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("% Auto-generated by ghandali_verify_index2.py\n")
        fh.write("\\begin{tabular}{llrrr}\n")
        fh.write("\\hline\n")
        fh.write("Category & Statement & Tested & Matches & Mismatches\\\\\n")
        fh.write("\\hline\n")
        for r in rows:
            fh.write(
                f"{latex_escape(r['category'])} & "
                f"{latex_escape(r['statement'])} & "
                f"{r['tested']} & {r['matches']} & {r['mismatches']}\\\\\n"
            )
        fh.write("\\hline\n")
        fh.write("\\end{tabular}\n")


# ============================================================
# XLSX writer, standard library only
# ============================================================

def _xml_escape(value):
    s = "" if value is None else str(value)
    return (
        s.replace("&", "&amp;")
         .replace("<", "&lt;")
         .replace(">", "&gt;")
         .replace('"', "&quot;")
    )


def _xlsx_col_name(idx):
    idx = int(idx)
    name = ""
    while True:
        idx, rem = divmod(idx, 26)
        name = chr(ord("A") + rem) + name
        if idx == 0:
            return name
        idx -= 1


def _xlsx_cell_xml(value, row_idx, col_idx):
    ref = f"{_xlsx_col_name(col_idx)}{row_idx}"

    if value is None:
        return f'<c r="{ref}"/>'

    if isinstance(value, bool):
        return f'<c r="{ref}"><v>{1 if value else 0}</v></c>'

    if isinstance(value, int):
        return f'<c r="{ref}"><v>{value}</v></c>'

    if isinstance(value, float) and math.isfinite(value):
        return f'<c r="{ref}"><v>{value:.17g}</v></c>'

    text = _xml_escape(value)
    return f'<c r="{ref}" t="inlineStr"><is><t>{text}</t></is></c>'


def _xlsx_sheet_xml(rows):
    out = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
        '<sheetData>',
    ]

    for r_idx, row in enumerate(rows, start=1):
        if row is None:
            row = []
        out.append(f'<row r="{r_idx}">')
        for c_idx, value in enumerate(row):
            out.append(_xlsx_cell_xml(value, r_idx, c_idx))
        out.append('</row>')

    out.extend([
        '</sheetData>',
        '</worksheet>',
    ])
    return "\n".join(out)


def _xlsx_safe_sheet_name(name, used):
    bad = set('[]:*?/\\')
    s = "".join("_" if ch in bad else ch for ch in str(name))
    s = s.strip() or "Sheet"
    s = s[:31]

    base = s
    k = 1
    while s in used:
        suffix = f"_{k}"
        s = (base[:31 - len(suffix)] + suffix)[:31]
        k += 1

    used.add(s)
    return s


def write_xlsx_workbook(path, sheets):
    path = Path(path)
    used = set()
    normalized = []

    for name, rows in sheets:
        safe = _xlsx_safe_sheet_name(name, used)
        normalized.append((safe, rows))

    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    ]

    for i in range(1, len(normalized) + 1):
        content_types.append(
            f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
            f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )

    content_types.append('</Types>')

    root_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1"
    Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument"
    Target="xl/workbook.xml"/>
</Relationships>
"""

    workbook_lines = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
        '<sheets>',
    ]

    for i, (name, _rows) in enumerate(normalized, start=1):
        workbook_lines.append(
            f'<sheet name="{_xml_escape(name)}" sheetId="{i}" r:id="rId{i}"/>'
        )

    workbook_lines.extend([
        '</sheets>',
        '</workbook>',
    ])

    workbook_rels = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    ]

    for i in range(1, len(normalized) + 1):
        workbook_rels.append(
            f'<Relationship Id="rId{i}" '
            f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{i}.xml"/>'
        )

    workbook_rels.append(
        f'<Relationship Id="rId{len(normalized) + 1}" '
        f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        f'Target="styles.xml"/>'
    )
    workbook_rels.append('</Relationships>')

    styles = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
  <fills count="1"><fill><patternFill patternType="none"/></fill></fills>
  <borders count="1"><border/></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>
  <cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>
"""

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", "\n".join(content_types))
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", "\n".join(workbook_lines))
        z.writestr("xl/_rels/workbook.xml.rels", "\n".join(workbook_rels))
        z.writestr("xl/styles.xml", styles)

        for i, (_name, rows) in enumerate(normalized, start=1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", _xlsx_sheet_xml(rows))


def _read_csv_rows(path):
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        for row in csv.reader(fh):
            rows.append(row)
    return rows


def _flatten_json(prefix, obj, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            _flatten_json(key, v, out)
    elif isinstance(obj, list):
        out.append([prefix, json.dumps(obj, ensure_ascii=False)])
    else:
        out.append([prefix, obj])


def _checks_to_rows(checks):
    rows = [["name", "tested", "matches", "mismatch", "first_mismatches", "value", "claim_status"]]
    for name, c in checks.items():
        if isinstance(c, dict):
            rows.append([
                name,
                c.get("tested", ""),
                c.get("matches", ""),
                c.get("mismatch", ""),
                json.dumps(c.get("first_mismatches", []), ensure_ascii=False),
                c.get("value", ""),
                c.get("claim_status", ""),
            ])
        else:
            rows.append([name, "", "", "", "", c, ""])
    return rows


def build_xlsx_from_existing_outputs(out_dir, xlsx_path):
    out_dir = Path(out_dir)
    xlsx_path = Path(xlsx_path)

    summary_path = out_dir / "verification-summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"missing {summary_path}")

    with open(summary_path, "r", encoding="utf-8") as fh:
        summary = json.load(fh)

    required = (
        "metadata",
        "theorem_checks",
        "exploratory_checks",
        "paper_table",
        "verdict",
    )
    missing = [k for k in required if k not in summary]
    if missing:
        raise KeyError(
            f"verification-summary.json missing keys {missing}; "
            f"found {sorted(summary)}"
        )

    metadata = summary.get("metadata", {})
    density = summary.get("density_context_not_theorem", {})
    theorem_checks = summary.get("theorem_checks", {})
    exploratory_checks = summary.get("exploratory_checks", {})
    paper_table = summary.get("paper_table", [])

    readme_rows = [
        ["key", "value"],
        ["title", metadata.get("title", "")],
        ["author", metadata.get("author", "")],
        ["script", metadata.get("script", "")],
        ["version", metadata.get("version", "")],
        ["generated_at", metadata.get("generated_at", "")],
        ["limit", metadata.get("limit", "")],
        ["sum_limit", metadata.get("sum_limit", "")],
        ["jac_limit", metadata.get("jac_limit", "")],
        ["workers", metadata.get("workers", "")],
        ["verdict", summary.get("verdict", "")],
        ["theorem_mismatch_total", summary.get("theorem_mismatch_total", "")],
        ["exploratory_mismatch_total", summary.get("exploratory_mismatch_total", "")],
        ["runtime_seconds", summary.get("runtime_seconds", "")],
        ["note", "Generated from verifier JSON/CSV outputs; theorem and exploratory checks are separated."],
    ]

    paper_rows = [["category", "name", "statement", "tested", "matches", "mismatches", "status"]]
    for r in paper_table:
        paper_rows.append([
            r.get("category", ""),
            r.get("name", ""),
            r.get("statement", ""),
            r.get("tested", ""),
            r.get("matches", ""),
            r.get("mismatches", ""),
            r.get("status", ""),
        ])

    density_rows = [["key", "value"]]
    _flatten_json("", density, density_rows)

    sheets = [
        ("README", readme_rows),
        ("paper_table", paper_rows),
        ("theorem_checks", _checks_to_rows(theorem_checks)),
        ("exploratory_checks", _checks_to_rows(exploratory_checks)),
        ("density_context", density_rows),
    ]

    numerical_csv = out_dir / "numerical-table.csv"
    if numerical_csv.exists():
        sheets.append(("numerical_table_csv", _read_csv_rows(numerical_csv)))

    records_csv = out_dir / "ghandali_index2_records.csv"
    if records_csv.exists():
        sheets.append(("records", _read_csv_rows(records_csv)))

    rank_csv = out_dir / "rank-index-all-primes.csv"
    if rank_csv.exists():
        sheets.append(("rank_inventory", _read_csv_rows(rank_csv)))

    log_path = out_dir / "verification-run.log"
    if log_path.exists():
        with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
            log_rows = [["line_no", "text"]]
            for i, line in enumerate(fh, start=1):
                log_rows.append([i, line.rstrip("\n")])
        sheets.append(("run_log", log_rows))

    write_xlsx_workbook(xlsx_path, sheets)


# ============================================================
# CLI
# ============================================================

def parse_index_set(s):
    if s is None or str(s).strip() == "":
        return set()
    out = set()
    for part in str(s).split(","):
        part = part.strip()
        if part:
            try:
                out.add(int(part))
            except ValueError:
                raise SystemExit(
                    f"--sum-indices: invalid entry {part!r}; expected comma-separated integers"
                )
    return out


def parse_args(argv):
    ap = argparse.ArgumentParser(
        description="Verifier for Fibonacci character sums at rank index two."
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    ap.add_argument("--interactive", action="store_true", help="force interactive setup UI")
    ap.add_argument("--limit", type=int, default=100000, help="prime limit N")
    ap.add_argument("--sum-limit", type=int, default=None, help="direct S(p) limit")
    ap.add_argument("--jac-limit", type=int, default=200000, help="direct Jacobsthal limit")
    ap.add_argument("--sum-indices", type=str, default="2,4", help="indices for direct S checks")
    ap.add_argument("--out", type=str, default=None,
                    help=f"output directory (default: {OFFICIAL_OUT}; with --safe: {SMOKE_OUT})")
    ap.add_argument("--workers", type=int, default=0, help="rank workers; 0 means auto CPU-2")
    ap.add_argument("--rank-chunk", type=int, default=5000, help="rank checkpoint chunk size")
    ap.add_argument("--checkpoint-every", type=int, default=25, help="checkpoint interval in chunks/records")
    ap.add_argument("--block", type=int, default=131072, help="Fibonacci block size")
    ap.add_argument("--jac-chunk", type=int, default=1000000, help="Jacobsthal chunk size")
    ap.add_argument("--no-sums", action="store_true", help="do not compute direct S(p)")
    ap.add_argument("--no-jac", action="store_true", help="do not compute Jacobsthal sums")
    ap.add_argument("--no-progress", action="store_true", help="disable one-line progress bars")
    ap.add_argument("--no-xlsx", action="store_true", help="do not write verification-workbook.xlsx")
    ap.add_argument("--xlsx-only", action="store_true", help="build verification-workbook.xlsx from existing outputs and exit")
    ap.add_argument("--fresh", action="store_true", help="discard checkpoint and start fresh")
    ap.add_argument("--resume", action="store_true", help="resume without asking if checkpoint exists")
    ap.add_argument("--safe", action="store_true", help=f"tiny safe smoke test (default output {SMOKE_OUT})")
    ap.add_argument("--force-official", action="store_true",
                    help="allow writing to the official directory for N < 2,000,000 or when it already holds a finished run")
    return ap.parse_args(argv)


# ============================================================
# Main
# ============================================================

def _same_dir(a, b):
    return Path(a).resolve() == Path(b).resolve()


def resolve_and_guard_output(args):
    """Choose the output directory and protect the official evidence directory."""
    if args.out is None:
        args.out = SMOKE_OUT if args.safe else OFFICIAL_OUT
    is_official = _same_dir(args.out, OFFICIAL_OUT)
    if args.safe and is_official:
        raise SystemExit(
            "ERROR: --safe must never write to the official evidence directory "
            f"({OFFICIAL_OUT}). Omit --out or use --out {SMOKE_OUT}."
        )
    if is_official and not args.xlsx_only and not args.force_official:
        if int(args.limit) < OFFICIAL_MIN_LIMIT:
            raise SystemExit(
                f"ERROR: the official directory {OFFICIAL_OUT} is reserved for the full run "
                f"(N >= {OFFICIAL_MIN_LIMIT}); got --limit {args.limit}. "
                "Use another --out, or --force-official."
            )
        finished = Path(args.out) / "verification-summary.json"
        if finished.exists():
            raise SystemExit(
                f"ERROR: {finished} already exists: the official directory holds a finished run. "
                "Move or delete the whole directory before a new official run (or use --force-official)."
            )
    return args


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    args = resolve_and_guard_output(args)
    
    if args.xlsx_only:
        out_dir = Path(args.out)
        xlsx_path = out_dir / "verification-workbook.xlsx"
        build_xlsx_from_existing_outputs(out_dir, xlsx_path)
        print(f"wrote {xlsx_path}  {xlsx_path.stat().st_size} bytes  sha-256 {sha256_file(xlsx_path)}")
        return 0

    if args.safe:
        args.limit = min(args.limit, 10000)
        args.sum_limit = 10000
        args.jac_limit = 10000
        args.workers = 1
        args.block = min(args.block, 32768)
        args.jac_chunk = min(args.jac_chunk, 100000)
        args.rank_chunk = min(args.rank_chunk, 1000)

    if args.sum_limit is None:
        args.sum_limit = args.limit

    if args.no_sums:
        args.sum_limit = 0
    if args.no_jac:
        args.jac_limit = 0

    args.limit = int(args.limit)
    args.sum_limit = int(args.sum_limit)
    args.jac_limit = int(args.jac_limit)
    args.rank_chunk = max(1, int(args.rank_chunk))
    args.checkpoint_every = max(1, int(args.checkpoint_every))

    if args.workers <= 0:
        args.workers = max(1, (cpu_count() or 3) - 2)

    sum_indices = parse_index_set(args.sum_indices)
    if args.no_sums:
        sum_indices = set()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    meta = make_checkpoint_meta(args, sum_indices)

    checkpoint_payload = None
    checkpoint_status = "not-requested"

    if args.fresh:
        clear_checkpoint(out_dir)
        checkpoint_status = "fresh-requested"
    else:
        loaded, status = load_checkpoint(out_dir, meta)
        checkpoint_status = status
        if loaded is not None:
            if args.resume:
                checkpoint_payload = loaded
            elif sys.stdin.isatty():
                if ask_resume_or_fresh(loaded):
                    checkpoint_payload = loaded
                else:
                    clear_checkpoint(out_dir)
                    checkpoint_status = "discarded-by-user"
            else:
                checkpoint_payload = loaded

    arch = archive_existing_outputs(out_dir)
    log = Logger(out_dir / "verification-run.log")

    t0 = time.time()

    log(f"{TITLE}")
    log(f"author: {AUTHOR}")
    log(f"script: ghandali_verify_index2.py {VERSION}")
    log(f"python: {platform.python_version()}  numpy: {np.__version__}")
    log(f"platform: {platform.platform()}")
    try:
        _sha = hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest()
    except OSError:
        _sha = "unavailable"
    log(f"script_sha256: {_sha}")
    log(f"argv: {sys.argv}")
    log(f"started: {_dt.datetime.now().isoformat(timespec='seconds')}")
    log(f"limit={args.limit}  sum_limit={args.sum_limit}  jac_limit={args.jac_limit}")
    log(f"workers={args.workers}  rank_chunk={args.rank_chunk}  sum_indices={args.sum_indices}")
    log(f"checkpoint_status={checkpoint_status}")
    if checkpoint_payload is not None:
        log("checkpoint loaded: resuming previous compatible run")
    if arch:
        log(f"archived previous outputs to: {arch}")

    progress_enabled = (not args.no_progress)

    # ---------- prime preparation ----------
    log("building prime sieve...")
    all_primes = primes_upto(args.limit)
    primes = [int(x) for x in all_primes if int(x) >= 3 and int(x) != 5]
    small_primes = [int(x) for x in primes_upto(math.isqrt(args.limit + 1) + 2)]
    log(f"primes for rank inventory: {len(primes)}")

    # ---------- restore checkpoint state ----------
    if checkpoint_payload is not None:
        rank_parts = {
            int(k): v for k, v in checkpoint_payload.get("rank_parts", {}).items()
        }
        recs = list(checkpoint_payload.get("recs", []))
    else:
        rank_parts = {}
        recs = []

    # ---------- rank inventory with checkpoint ----------
    log("computing ranks and indices...")

    chunks = []
    for idx, start in enumerate(range(0, len(primes), args.rank_chunk)):
        chunks.append((idx, primes[start:start + args.rank_chunk]))

    done_rank = set(int(k) for k in rank_parts.keys())
    pending_chunks = [
        (idx, ch) for idx, ch in chunks
        if idx not in done_rank
    ]

    rank_bar = ProgressBar(
        "rank chunks",
        len(chunks),
        enabled=progress_enabled,
    )
    rank_bar.update(len(done_rank), note="resume" if done_rank else "", force=True)

    tasks_done_since_ckpt = 0

    def collect_rank_result(idx, rows):
        nonlocal tasks_done_since_ckpt
        rank_parts[int(idx)] = np.array(rows, dtype=np.int64)
        tasks_done_since_ckpt += 1
        done_now = len(rank_parts)
        p_note = ""
        if rows:
            p_note = f"p≤{rows[-1][0]}"
        rank_bar.update(done_now, note=p_note)

        if tasks_done_since_ckpt >= args.checkpoint_every:
            save_checkpoint(out_dir, meta, rank_parts, recs, phase="rank")
            tasks_done_since_ckpt = 0

    try:
        if pending_chunks:
            if args.workers == 1:
                for idx, ch in pending_chunks:
                    ridx, rows = _rank_worker((idx, ch, small_primes))
                    collect_rank_result(ridx, rows)
            else:
                payloads = [(idx, ch, small_primes) for idx, ch in pending_chunks]
                with Pool(processes=args.workers) as pool:
                    for ridx, rows in pool.imap_unordered(_rank_worker, payloads):
                        collect_rank_result(ridx, rows)

        save_checkpoint(out_dir, meta, rank_parts, recs, phase="rank-complete")
        rank_bar.finish("rank complete")

    except KeyboardInterrupt:
        rank_bar.finish("interrupted")
        save_checkpoint(out_dir, meta, rank_parts, recs, phase="rank-interrupted")
        log("Interrupted during rank stage. Checkpoint saved.")
        log.close()
        return 130

    if rank_parts:
        R = np.concatenate([rank_parts[i] for i in sorted(rank_parts)], axis=0)
        R = R[np.argsort(R[:, 0])]
    else:
        R = np.zeros((0, 5), dtype=np.int64)

    log(f"rank rows: {len(R)}")

    # ---------- direct S records with checkpoint ----------
    recs_by_p = {}
    for r in recs:
        try:
            recs_by_p[int(r["p"])] = r
        except Exception:
            pass

    if not args.no_sums and args.sum_limit > 0 and sum_indices:
        selected = [
            row for row in R
            if int(row[0]) <= args.sum_limit
            and int(row[3]) in sum_indices
            and int(row[0]) >= 7
        ]

        selected_ps = set(int(row[0]) for row in selected)
        recs_by_p = {
            p: r for p, r in recs_by_p.items()
            if p in selected_ps
        }

        pending_rows = [
            row for row in selected
            if int(row[0]) not in recs_by_p
        ]

        log(f"direct S records selected: {len(selected)}")
        if recs_by_p:
            log(f"direct S records already in checkpoint: {len(recs_by_p)}")
        log("computing direct S(p) records...")

        s_bar = ProgressBar(
            "S records",
            len(selected),
            enabled=progress_enabled,
        )
        s_bar.update(len(recs_by_p), note="resume" if recs_by_p else "", force=True)

        done_since_ckpt = 0

        try:
            for row in pending_rows:
                p_now = int(row[0])
                rec = build_sum_record(row, small_primes, args)
                recs_by_p[p_now] = rec
                done_since_ckpt += 1

                s_bar.update(len(recs_by_p), note=f"p={p_now}")

                if done_since_ckpt >= args.checkpoint_every:
                    recs = [recs_by_p[p] for p in sorted(recs_by_p)]
                    save_checkpoint(out_dir, meta, rank_parts, recs, phase="sums")
                    done_since_ckpt = 0

            recs = [recs_by_p[p] for p in sorted(recs_by_p)]
            save_checkpoint(out_dir, meta, rank_parts, recs, phase="sums-complete")
            s_bar.finish("S complete")

        except KeyboardInterrupt:
            s_bar.finish("interrupted")
            recs = [recs_by_p[p] for p in sorted(recs_by_p)]
            save_checkpoint(out_dir, meta, rank_parts, recs, phase="sums-interrupted")
            log("Interrupted during direct S stage. Checkpoint saved.")
            log.close()
            return 130

    else:
        log("direct S records skipped.")
        recs = [recs_by_p[p] for p in sorted(recs_by_p)]

    # ---------- aggregate ----------
    recs.sort(key=lambda r: r["p"])
    density, theorem, exploratory, max_ratio = aggregate_checks(R, recs)
    table = build_paper_table(theorem, exploratory)

    theorem_mismatch_total = sum(
        int(v.get("mismatch", 0))
        for v in theorem.values()
        if isinstance(v, dict) and "mismatch" in v
    )
    exploratory_mismatch_total = sum(
        int(v.get("mismatch", 0))
        for v in exploratory.values()
        if isinstance(v, dict) and "mismatch" in v
    )

    runtime = time.time() - t0
    verdict = "PASS" if theorem_mismatch_total == 0 else "FAIL"

    summary = {
        "metadata": {
            "title": TITLE,
            "author": AUTHOR,
            "script": "ghandali_verify_index2.py",
            "version": VERSION,
            "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
            "limit": args.limit,
            "sum_limit": args.sum_limit,
            "jac_limit": args.jac_limit,
            "workers": args.workers,
            "rank_chunk": args.rank_chunk,
            "checkpoint_every": args.checkpoint_every,
            "sum_indices": sorted(sum_indices),
            "sums_enabled": not args.no_sums,
            "jacobsthal_enabled": not args.no_jac,
            "xlsx_enabled": not args.no_xlsx,
        },
        "verdict": verdict,
        "theorem_mismatch_total": theorem_mismatch_total,
        "exploratory_mismatch_total": exploratory_mismatch_total,
        "runtime_seconds": runtime,
        "density_context_not_theorem": density,
        "theorem_checks": theorem,
        "exploratory_checks": exploratory,
        "paper_table": table,
    }

    # ---------- write outputs ----------
    rank_path = out_dir / "rank-index-all-primes.csv"
    records_path = out_dir / "ghandali_index2_records.csv"
    table_csv = out_dir / "numerical-table.csv"
    table_tex = out_dir / "numerical-table.tex"
    summary_path = out_dir / "verification-summary.json"
    xlsx_path = out_dir / "verification-workbook.xlsx"

    write_rank_csv(R, rank_path)
    write_records_csv(recs, records_path)
    write_table_csv(table, table_csv)
    write_table_tex(table, table_tex)

    with open(summary_path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
        fh.write("\n")

    # Successful completion: remove checkpoint.
    clear_checkpoint(out_dir)

    log("")
    log("THEOREM CHECKS")
    for name, c in theorem.items():
        if isinstance(c, dict) and "mismatch" in c:
            log(
                f"  {name:65s} tested={c['tested']:7d} "
                f"mismatch={c['mismatch']:4d}"
            )

    log("")
    log("EXPLORATORY CHECKS")
    for name, c in exploratory.items():
        if isinstance(c, dict) and "mismatch" in c:
            log(
                f"  {name:65s} tested={c['tested']:7d} "
                f"mismatch={c['mismatch']:4d}"
            )
        elif isinstance(c, dict) and "value" in c:
            log(f"  {name:65s} value={c['value']}")

    log("")
    if density:
        log(
            "density context: "
            f"index2 split/A={density.get('index2_split_over_A', 0.0):.6f} "
            f"(context 9/38={9/38:.6f}), "
            f"index2 inert/A={density.get('index2_inert_over_A', 0.0):.6f} "
            f"(context 10/19={10/19:.6f}, Dirichlet-conditional)"
        )

    log("")
    for pth in (rank_path, records_path, table_csv, table_tex, summary_path):
        log(
            f"wrote {pth.name:36s} "
            f"{pth.stat().st_size:10d} bytes  sha-256 {sha256_file(pth)}"
        )

    log("")
    log("checkpoint cleared after successful completion")
    log(f"VERDICT: {verdict}")
    log(f"theorem mismatches: {theorem_mismatch_total}")
    log(f"exploratory mismatches: {exploratory_mismatch_total}")
    log(f"runtime seconds: {runtime:.2f}")
    log.close()

    if not args.no_xlsx:
        try:
            build_xlsx_from_existing_outputs(out_dir, xlsx_path)
            print(
                f"wrote {xlsx_path.name}  {xlsx_path.stat().st_size} bytes  "
                f"sha-256 {sha256_file(xlsx_path)}"
            )
        except Exception as exc:
            print(f"WARNING: xlsx build failed (non-fatal): {exc}")

    print(f"wrote {log.path.name}  sha-256 {sha256_file(log.path)}")

    try:
        shutil.copy2(Path(__file__).resolve(), out_dir / "ghandali_verify_index2.py")
    except OSError as exc:
        print(f"WARNING: could not copy verifier into output directory: {exc}")

    try:
        man = write_manifest(out_dir)
        print(f"wrote manifest {man}  sha-256 {sha256_file(man)}")
    except OSError as exc:
        print(f"WARNING: could not write sha256sums.txt: {exc}")

    return 0 if theorem_mismatch_total == 0 else 1


if __name__ == "__main__":
    freeze_support()
    argv_rest = [a for a in sys.argv[1:] if a != "--interactive"]
    want_ui = (
        (not argv_rest)
        or ("--interactive" in sys.argv[1:])
    )
    if want_ui and sys.stdin.isatty():
        sys.argv = [sys.argv[0]] + interactive_setup_argv()
    elif want_ui and not sys.stdin.isatty():
        sys.argv = [sys.argv[0]] + argv_rest
    else:
        sys.argv = [sys.argv[0]] + argv_rest
    raise SystemExit(main())