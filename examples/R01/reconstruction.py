"""Exact Figure 3 reconstruction from anonymised R01 workspace values.

Importing this module computes the two arrays without writing output files.
The fixed mapping is the documented post hoc mapping, not a new diagnostic.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
HERE = Path(__file__).resolve().parent
CH, PIX, V = 'Fe56', 15.0, 375.0
BASE = float(2**30)  # Synthetic arithmetic origin, not a real timestamp.
rep = json.loads((HERE / "reference_mapping.json").read_text(encoding="utf-8"))
blk = pd.read_csv(HERE / 'reference_blocks.csv', float_precision='round_trip')
for _k in ('t_mid', 'observed_start_s', 'observed_end_s'):
    blk[_k] += BASE
assert rep["status"] == "supported" and rep["channel"] == CH
ALIGN = rep["alignment"]["offset_s"]; R_IT = rep["rate_ppm"] * 1e-6; RES = rep["closure"]["residual_fit"]["intercept_s"]
a = np.loadtxt(HERE / 'R01_signals.csv', delimiter=',', skiprows=1)
t = a[:, 0] + BASE
y = a[:, 2]
rows = pd.read_csv(HERE / 'R01_workspace_intervals.csv', float_precision='round_trip').to_dict(orient='records')
L = [dict(**{'Start Time': d['start_s'] + BASE, 'End Time': d['end_s'] + BASE},
          StartX=d['start_x_um'], EndX=d['end_x_um'], StartY=d['start_y_um'], EndY=d['end_y_um']) for d in rows]
assert len(L) == rep["n_lines"] == 904
ls = np.array([d["Start Time"] for d in L]); le = np.array([d["End Time"] for d in L])
sx = np.array([d["StartX"] for d in L]); ex = np.array([d["EndX"] for d in L]); ys = np.array([d["StartY"] for d in L])
T0 = ls[0]; mid = (ls + le) / 2
model = lambda tt: ALIGN + RES + R_IT * (tt - T0)
tau_line = model(mid); tau_c = float(np.median(tau_line))
disp_um = (tau_line - tau_c) * V
el_h = (mid - T0) / 3600.0
tau_blk = ALIGN + blk.lag_s.to_numpy(); tmid_blk = blk.t_mid.to_numpy()
blk_h = (tmid_blk - T0) / 3600.0
retained = np.array([b in set(rep["closure"].get("retained_blocks") or []) for b in blk.block])

x0 = float(min(sx.min(), ex.min())); ncol = int(np.ceil((max(sx.max(), ex.max()) - x0) / PIX)) + 1
uy = np.unique(np.round(ys, 1)); rowof = {val: i for i, val in enumerate(uy)}


def raster(lag_s, lag_e):
    """Reconstruct with lag evaluated at each line endpoint; no intensity resampling before binning."""
    s, e = ls + lag_s, le + lag_e
    idx = np.searchsorted(s, t, side="right") - 1
    ok = idx >= 0
    ic = np.clip(idx, 0, len(L) - 1)
    ok &= t <= e[ic]
    li, tt, yy = ic[ok], t[ok], y[ok]
    xp = sx[li] + (tt - s[li]) / (e[li] - s[li]) * (ex[li] - sx[li])
    col = np.clip(np.floor((xp - x0) / PIX).astype(int), 0, ncol - 1)
    row = np.array([rowof[val] for val in np.round(ys[li], 1)])
    acc = np.zeros((len(uy), ncol)); cnt = np.zeros_like(acc)
    np.add.at(acc, (row, col), yy); np.add.at(cnt, (row, col), 1)
    with np.errstate(invalid="ignore"):
        img = acc / cnt
    img[cnt == 0] = np.nan
    for r in range(img.shape[0]):
        vv = np.flatnonzero(np.isfinite(img[r]))
        if len(vv) >= 2:
            ins = np.arange(vv.min(), vv.max() + 1)
            img[r, ins] = np.interp(ins, vv, img[r, vv])
    return img


q = lambda img: np.where(np.isfinite(img) & (img > 0), np.log10(np.clip(img, 1, None)), np.nan)
qc = q(raster(np.full(len(L), tau_c), np.full(len(L), tau_c))); ql = q(raster(model(ls), model(le)))
ql_mid = q(raster(tau_line, tau_line))                     # Midpoint-shift approximation, retained only for comparison.
_both = np.isfinite(ql) & np.isfinite(ql_mid)
EXACT_VS_MID = dict(cells_compared=int(_both.sum()), cells_different=int(np.sum(_both & (np.abs(ql - ql_mid) > 1e-12))),
                    finite_mask_differs=int(np.sum(np.isfinite(ql) != np.isfinite(ql_mid))),
                    max_abs_log10_diff=float(np.nanmax(np.abs(ql - ql_mid)[_both])) if _both.any() else 0.0,
                    max_line_end_shift_um=float(np.max(R_IT * (le - ls) / 2) * V))   # Maximum endpoint-versus-midpoint shift, in micrometres.
both = np.isfinite(qc) & np.isfinite(ql)
lo, hi = np.nanpercentile(np.r_[qc[both], ql[both]], [2, 98])
ext = [0, ncol * PIX / 1000.0, (uy.max() - uy.min()) / 1000.0, 0]
yrow = (uy - uy.min()) / 1000.0
row_h = np.interp(uy, np.sort(ys), el_h[np.argsort(ys)])

WIN = {"B": dict(x=(9.0, 15.0), y=(1.2, 3.6), name="Early"), "C": dict(x=(0.5, 6.5), y=(9.6, 12.0), name="Late")}
for k, w in WIN.items():
    m = (yrow >= w["y"][0]) & (yrow <= w["y"][1])
    w["h"] = (float(row_h[m].min()), float(row_h[m].max()))
    lm = (el_h >= w["h"][0]) & (el_h <= w["h"][1])
    w["shift_mm"] = float(np.mean(disp_um[lm])) / 1000
    bm = retained & (blk_h >= w["h"][0]) & (blk_h <= w["h"][1])
    w["blocks"] = [int(b) for b in blk.block[bm]]
    w["rms_const_mm"] = float(np.sqrt(np.mean((tau_blk[bm] - tau_c) ** 2)) * V / 1000)
    w["rms_corr_mm"] = float(np.sqrt(np.mean((tau_blk[bm] - model(tmid_blk[bm])) ** 2)) * V / 1000)

