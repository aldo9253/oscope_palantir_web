#!/usr/bin/env python3
"""
Quick helper to convert a folder of .trc files into CSV files grouped
by event and channel for use with the web UI.

Usage:
  python trc_to_csv.py /path/to/trc_folder /path/to/output_csv

Requires `readTrcDoner.py` available on the Python path.
"""
import sys, os, re, csv
from pathlib import Path

try:
    from readTrcDoner import Trc
except Exception as e:
    print("ERROR: Could not import readTrcDoner.Trc:", e)
    sys.exit(1)

def main(inp, outdir):
    trc = Trc()
    os.makedirs(outdir, exist_ok=True)
    files = [f for f in os.listdir(inp) if f.lower().endswith('.trc')]
    pat = re.compile(r"C(\d+).*-([0-9]+)\.trc$", re.I)
    if not files:
        print("No .trc files found in", inp)
        return
    for f in files:
        m = pat.search(f)
        if not m:
            continue
        ch = int(m.group(1))
        evt = m.group(2)
        path = os.path.join(inp, f)
        try:
            t, y, meta = trc.open(path)
        except Exception as e:
            print("Failed to read", f, e)
            continue
        # Detrend by removing mean, like the desktop app
        ym = float(sum(y)/len(y)) if len(y)>0 else 0.0
        y = [float(v)-ym for v in y]
        out = Path(outdir) / f"C{ch}-{evt}.csv"
        with open(out, 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['time','voltage'])
            for tt,vv in zip(t,y):
                w.writerow([tt, vv])
        print("Wrote", out)

if __name__ == '__main__':
    if len(sys.argv)<3:
        print("Usage: python trc_to_csv.py INPUT_TRC_FOLDER OUTPUT_CSV_FOLDER")
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])

