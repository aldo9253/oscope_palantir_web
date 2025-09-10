Oscilloscope Palantir (Web)
===========================

This is a browser-based reimplementation (in progress) of `oscope_palantir.py`.
It aims to provide a similar UI and workflow using HTML/JavaScript.

Status
------
- Loads local waveform files (CSV) via the browser (no server required).
- Loads LeCroy `.trc` files directly in the browser (experimental) — no CSV expansion required.
- Groups files into events by filename and channel.
- Plots selected event channels on an interactive canvas.
- Savitzky–Golay filtering (SG) for smoothing.
- Dynamic fits: QD3Fit, QDMFit, CSA_pulse, skew_gaussian, gaussian.
- Batch Run across an event range for the selected fit/time window.
- Feature Scan to flag events with |signal| > threshold × std.
- Results are kept in-memory and can be exported as JSON; CSV fallback export is available.
- Import Results (JSON) merges saved fits back into the session (dataset-aware).
- low_pass_max measurement (SG+peak) with invert handling.
- Derived metrics for QD3/QDM (charge/mass/radius) shown in Fit Info.

Limitations (for now)
---------------------
- `.trc` support is based on the LeCroy LECROY_2_3 WAVEDESC template used by
  `readTrcDoner.py`. Other templates may not parse correctly yet.
- CSV input is still supported as an alternative for testing.

Usage
-----
1) Open `index.html` in a modern browser (Chrome/Edge/Firefox).
2) Click “Select Folder” and choose your CSV folder (hold Shift/Ctrl to select multiple
   files, or choose a directory if your browser supports `webkitdirectory`).
3) Choose an event from the dropdown to plot.
4) Use SG filter, run dynamic fits, batch run, and feature scan similarly to the desktop app.

File naming convention
----------------------
Files should be named like `C<channel>-<event>.<ext>` (e.g., `C1-00012.trc`,
`C2-00012.csv`). The app groups files with the same `<event>` number into one
event, per channel. A more flexible pattern `C<ch>...-<event>.trc` is also
recognized for `.trc` files.

Helper: TRC to CSV (optional)
------------------------------
Use `trc_to_csv.py` to convert `.trc` files into CSV pairs per event:

    python trc_to_csv.py /path/to/trc_folder /path/to/output_csv

This requires Python with `readTrcDoner.py` available in your PYTHONPATH.

Exporting results
-----------------
Click “Export Results” to download a JSON of fit results. This mirrors the in-memory
structure with keys like `evt_<id>_<fit>_ch<ch>`.

Planned
-------
- Direct `.trc` parsing in the browser (if feasible), or a lightweight local companion.
- HDF5 export/import via h5wasm (scaffolded; see below to enable).
- UI polish to fully match the desktop app.

HDF5 (optional)
---------------
To enable HDF5 import/export in the browser, place the h5wasm bundle under:

    oscope_html/vendor/h5wasm/hdf5_hl.js

and the associated WASM assets as required by that build. After reloading
`index.html`, the Export/Import HDF5 buttons will attempt to use the HDF5 engine.
If unavailable, the app offers a CSV fallback for export and a message for import.
