#!/usr/bin/env python3
"""h9_csv — add Hex9 addresses to a CSV of points, for dataset prep.

Reads a CSV that has longitude/latitude columns and writes the same CSV with an
appended `h9_uuid` column holding the full (reversible) Hex9 UUID per row — the
self-contained address you store, index, and JOIN against h9_grid downstream.
Optionally appends a canonical bin column (`--bin L`) and/or a human label
column (`--label L`). Streams in chunks, so it handles large files in bounded
memory, and batch-encodes (OpenMP) for speed — suitable for background jobs that
prepare datasets before they land in SQL.

Examples
--------
  # auto-detect lon/lat columns, write to stdout
  h9_csv.py points.csv > points_h9.csv

  # explicit columns, add an L12 canonical bin column too
  h9_csv.py --lon x --lat y --bin 12 in.csv out.csv

  # from a pipe
  cat points.csv | h9_csv.py - - > out.csv

The full UUID is the load-bearer: re-bin or label it at any layer later. Bins
and labels are derived/canonical (grid-matching); see the project docs.
"""
from __future__ import annotations

import argparse
import binascii
import csv
import os
import sys

# Locate the hex9_ext extension (built into ../build by default).
_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, os.pardir, "build"), os.path.join(_HERE, os.pardir)):
    if os.path.isdir(_cand):
        sys.path.insert(0, os.path.abspath(_cand))
try:
    import numpy as np
    import hex9_ext as h9
except ImportError as e:  # pragma: no cover - environment guidance
    sys.exit(f"h9_csv: cannot import hex9_ext / numpy ({e}).\n"
             f"  Build the extension first (cmake --build build) or set "
             f"PYTHONPATH to its directory.")

LON_NAMES = ("longitude", "lon", "lng", "long", "x")
LAT_NAMES = ("latitude", "lat", "y")


def _find_col(header, wanted, names, role):
    """Resolve a column name to an index, by explicit name or autodetection."""
    lower = [c.strip().lower() for c in header]
    if wanted is not None:
        if wanted in header:
            return header.index(wanted)
        if wanted.lower() in lower:
            return lower.index(wanted.lower())
        sys.exit(f"h9_csv: {role} column {wanted!r} not found in header {header}")
    for n in names:
        if n in lower:
            return lower.index(n)
    sys.exit(f"h9_csv: could not autodetect the {role} column (looked for "
             f"{', '.join(names)}); pass --{role} explicitly. Header: {header}")


def _uuids_from_bytes(arr) -> list[str]:
    """(m,16) uint8 -> list of hyphenated UUID strings (vectorised hexlify)."""
    hexs = binascii.hexlify(np.ascontiguousarray(arr).tobytes()).decode("ascii")
    out = []
    for i in range(0, len(hexs), 32):
        s = hexs[i:i + 32]
        out.append(f"{s[0:8]}-{s[8:12]}-{s[12:16]}-{s[16:20]}-{s[20:32]}")
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="h9_csv",
        description="Append Hex9 addresses (full UUID, optional bin/label) to a CSV.")
    p.add_argument("input", nargs="?", default="-",
                   help="input CSV path, or - for stdin (default)")
    p.add_argument("output", nargs="?", default="-",
                   help="output CSV path, or - for stdout (default)")
    p.add_argument("--lon", help="longitude column name (default: autodetect)")
    p.add_argument("--lat", help="latitude column name (default: autodetect)")
    p.add_argument("--out-col", default="h9_uuid",
                   help="name of the full-UUID column to add (default: h9_uuid)")
    p.add_argument("--bin", type=int, metavar="LAYER",
                   help="also add a canonical bin column h9_bin<LAYER> (0..29)")
    p.add_argument("--label", type=int, metavar="LAYER",
                   help="also add a label column h9_label<LAYER> (0..29)")
    p.add_argument("--delimiter", default=",", help="CSV delimiter (default: ,)")
    p.add_argument("--chunk", type=int, default=100_000,
                   help="rows per batch (default: 100000)")
    args = p.parse_args(argv)

    for L in (args.bin, args.label):
        if L is not None and not (0 <= L <= 29):
            sys.exit("h9_csv: --bin/--label layer must be 0..29")

    h9.warp_init()

    fin = sys.stdin if args.input == "-" else open(args.input, newline="")
    fout = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    try:
        reader = csv.reader(fin, delimiter=args.delimiter)
        try:
            header = next(reader)
        except StopIteration:
            sys.exit("h9_csv: empty input")
        lon_i = _find_col(header, args.lon, LON_NAMES, "lon")
        lat_i = _find_col(header, args.lat, LAT_NAMES, "lat")

        extra = [args.out_col]
        if args.bin is not None:
            extra.append(f"h9_bin{args.bin}")
        if args.label is not None:
            extra.append(f"h9_label{args.label}")
        writer = csv.writer(fout, delimiter=args.delimiter)
        writer.writerow(header + extra)

        total = skipped = 0

        def flush(rows):
            nonlocal skipped
            lons, lats, valid = [], [], []
            for j, row in enumerate(rows):
                try:
                    lons.append(float(row[lon_i]))
                    lats.append(float(row[lat_i]))
                    valid.append(j)
                except (ValueError, IndexError):
                    skipped += 1
            cols = [[""] * len(rows) for _ in extra]
            if valid:
                full = h9.encode(np.asarray(lons), np.asarray(lats))
                series = [_uuids_from_bytes(full)]                       # full UUID
                if args.bin is not None:
                    series.append(_uuids_from_bytes(h9.bin(full, args.bin)))
                if args.label is not None:
                    series.append([h9.label(row, args.label, False) for row in full])
                for ci, vals in enumerate(series):
                    for k, j in enumerate(valid):
                        cols[ci][j] = vals[k]
            for j, row in enumerate(rows):
                writer.writerow(row + [c[j] for c in cols])

        batch = []
        for row in reader:
            batch.append(row)
            if len(batch) >= args.chunk:
                total += len(batch)
                flush(batch)
                batch = []
        if batch:
            total += len(batch)
            flush(batch)
    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    sys.stderr.write(f"h9_csv: wrote {total} rows"
                     + (f", {skipped} skipped (unparseable lon/lat)" if skipped else "")
                     + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
