"""Wheel-tag honesty gate: an abi3-tagged wheel must contain an abi3 extension.

hex9 2.2.0–2.3.0 shipped cp312-abi3 wheels holding a version-specific
_core.cpython-312-*.so (nanobind silently drops STABLE_ABI when
Development.SABIModule isn't found) — such a wheel installs on every
CPython >= 3.12 and then fails to import on all but 3.12
(found by dggs_compare CI, ajfriend/dggs_compare#35).

Usage: python test/abi3_check.py <wheelhouse-dir>
Scans every *.whl; for wheels whose tag claims abi3, the bundled extension
must be version-agnostic (_core.abi3.so / _core.pyd), never cpython-3XX/cp3XX.
"""
import re
import sys
import zipfile
from pathlib import Path

EXT_RE = re.compile(r"\.(so|pyd|dylib)$")
VERSIONED_RE = re.compile(r"(cpython-3\d+|\.cp3\d+-)")


def check(wheel: Path) -> list[str]:
    errors = []
    with zipfile.ZipFile(wheel) as z:
        exts = [n for n in z.namelist() if EXT_RE.search(n)]
    if "abi3" in wheel.name:
        if not exts:
            errors.append(f"{wheel.name}: abi3 tag but no extension module inside")
        for n in exts:
            if VERSIONED_RE.search(n):
                errors.append(
                    f"{wheel.name}: abi3 tag but version-specific extension {n} "
                    "(STABLE_ABI silently dropped? see test/abi3_check.py)"
                )
    return errors


def main() -> int:
    wheelhouse = Path(sys.argv[1] if len(sys.argv) > 1 else "wheelhouse")
    wheels = sorted(wheelhouse.glob("*.whl"))
    if not wheels:
        print(f"abi3_check: no wheels found in {wheelhouse}", file=sys.stderr)
        return 1
    errors = [e for w in wheels for e in check(w)]
    for e in errors:
        print(f"abi3_check: FAIL {e}", file=sys.stderr)
    if not errors:
        print(f"abi3_check: OK ({len(wheels)} wheel(s), tags match contents)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
