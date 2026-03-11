Oscilloscope Viewer (Web)
=========================

Browser-based waveform viewer for local oscilloscope data (`.trc`, `.csv`, `.txt`).

What It Does
------------
- Loads local waveform files directly in the browser.
- Parses LeCroy `.trc` files in-browser (worker + fallback parser).
- Groups files by event using filename patterns like `C1-00012.trc`.
- Displays event waveforms per channel with pan/zoom and decimation.
- Supports SG-filter overlay for quick smoothing checks.
- Lets you show/hide channels with a channel visibility menu:
  - `Show All`
  - `Hide All`
  - Per-channel toggles

Notes
-----
- This web app is now viewing-focused.
- Fit routines, fit/batch-fit tools, and fit result export/import workflows were removed.

Usage
-----
1. Open `index.html` in a modern browser.
2. Click `Select Folder` and choose your data folder.
3. Pick an event from the dropdown.
4. Use:
   - `Tool: Zoom` to drag a zoom box
   - `Tool: Pan` to drag view (hold `Shift` for vertical pan)
   - mouse wheel to zoom (`Shift`/`Ctrl` for Y zoom)
   - double-click a plot to reset that channel view
5. Use `Channels` menu to hide/show channel plots.
6. Optional: run `SG Filter` on a selected channel.

File Naming
-----------
Supported grouping patterns:
- `C<channel>-<event>.<ext>`
- `C<channel>...-<event>.trc`

Examples:
- `C1-00012.trc`
- `C2-00012.csv`
