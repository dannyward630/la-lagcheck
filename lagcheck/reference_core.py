# Extracted unchanged from study estimator.
# Source SHA256: 438063d3fd0a8794264943ac254ab429a324884f68ff0b90b5ddeb57ecfe3e22
import numpy as np
import pandas as pd
MAX_SHIFT=25.0
N_BLOCKS=10
MIN_LINES_PER_BLOCK=12
PEAK_MIN=0.25
MARGIN_MIN=0.01

def on_sequences(t, sig, lines):
    """Binary laser-on from the log and signal-on from the mass spectrometer, on the MS grid."""
    log_on = np.zeros(len(t), dtype=np.float64)
    for d in lines:
        i0 = np.searchsorted(t, d["Start Time"])
        i1 = np.searchsorted(t, d["End Time"])
        if i1 > i0:
            log_on[i0:i1] = 1.0
    # Matched filter: the log side is binary (laser on / off), the signal side keeps its
    # full log-intensity waveform. Binary-vs-binary loses contrast when the duty cycle is
    # high, which is the case for tissue sections where only the short inter-line gap is
    # truly off.
    v = np.log10(np.clip(np.nan_to_num(np.asarray(sig, float), nan=1.0), 1, None))
    lo, hi = np.percentile(v, [5, 95])
    rng = max(hi - lo, 1e-6)
    sig_wave = np.clip((v - lo) / rng, -0.5, 1.5)
    return log_on, sig_wave

def best_shift(a, b, maxlag, dt):
    """Pearson correlation of a against b over integer shifts; returns (lag_s, peak, margin).

    a is the log pattern for one block, b the signal pattern over the same span widened by
    maxlag on each side. margin is the peak minus the best value more than 1 s away, so a
    peak that barely beats its surroundings can be rejected.
    """
    n = len(a)
    a = a - a.mean()
    na = np.sqrt((a * a).sum())
    if na <= 0:
        return np.nan, np.nan, np.nan
    # rolling mean and norm of every window of b of length n
    cs = np.concatenate([[0.0], np.cumsum(b)])
    cs2 = np.concatenate([[0.0], np.cumsum(b * b)])
    m = len(b) - n + 1
    if m < 3:
        return np.nan, np.nan, np.nan
    s1 = cs[n:] - cs[:-n]
    s2 = cs2[n:] - cs2[:-n]
    mean_b = s1 / n
    var_b = np.maximum(s2 - n * mean_b**2, 0.0)
    nb = np.sqrt(var_b)
    # dot products for every alignment, via FFT
    size = 1
    while size < len(b) + n:
        size *= 2
    fa = np.fft.rfft(a[::-1], size)
    fb = np.fft.rfft(b, size)
    conv = np.fft.irfft(fa * fb, size)
    dots = conv[n - 1: n - 1 + m]
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = dots / (na * nb)
    corr[~np.isfinite(corr)] = -np.inf
    k = int(np.argmax(corr))
    # parabolic sub-sample interpolation of the peak, so lags are not quantised to the
    # sweep grid (without it every lag is an integer multiple of dt and the residual
    # scatter sits at the dt/sqrt(12) quantisation floor)
    frac = 0.0
    if 0 < k < m - 1 and np.isfinite(corr[k - 1]) and np.isfinite(corr[k + 1]):
        denom = corr[k - 1] - 2.0 * corr[k] + corr[k + 1]
        if denom < 0:
            frac = float(np.clip(0.5 * (corr[k - 1] - corr[k + 1]) / denom, -0.5, 0.5))
    lag_samples = k + frac - round(maxlag / dt)
    far = np.abs(np.arange(m) - k) * dt > 1.0
    margin = float(corr[k] - corr[far].max()) if far.any() else np.nan
    return lag_samples * dt, float(corr[k]), margin

def analyse_loaded(t, sig, L, sample, n_blocks=None, max_shift=None,
                   peak_min=None, margin_min=None, reject_tol=3.0):
    """The estimator itself, on arrays already in memory.

    Rate-recovery and sensitivity checks use this same code path with modified
    log times or settings.
    """
    N_BLOCKS = globals()["N_BLOCKS"] if n_blocks is None else int(n_blocks)
    MAX_SHIFT = globals()["MAX_SHIFT"] if max_shift is None else float(max_shift)
    PEAK_MIN = globals()["PEAK_MIN"] if peak_min is None else float(peak_min)
    MARGIN_MIN = globals()["MARGIN_MIN"] if margin_min is None else float(margin_min)
    if len(L) < N_BLOCKS * MIN_LINES_PER_BLOCK:
        return None, None

    # Mean spacing, not median: the MS time grid is non-uniform with a quantised, skewed
    # spacing distribution (median exceeds mean by up to 6000 ppm on some runs), and the
    # index-shift-to-seconds conversion needs the average spacing.
    dt = float((t[-1] - t[0]) / (len(t) - 1))
    log_on, sig_on = on_sequences(t, sig, L)
    pad = int(round(MAX_SHIFT / dt))

    ls = np.array([d["Start Time"] for d in L])
    edges = np.linspace(0, len(L), N_BLOCKS + 1).astype(int)
    recs = []
    for bi in range(N_BLOCKS):
        seg = L[edges[bi]: edges[bi + 1]]
        if len(seg) < MIN_LINES_PER_BLOCK:
            continue
        i0 = np.searchsorted(t, seg[0]["Start Time"])
        i1 = np.searchsorted(t, seg[-1]["End Time"])
        if i1 - i0 < 200:
            continue
        a = log_on[i0:i1]
        j0, j1 = max(0, i0 - pad), min(len(t), i1 + pad)
        b = sig_on[j0:j1]
        if len(b) < len(a) + 3:
            continue
        # shift is measured relative to the window actually used
        lag, peak, margin = best_shift(a, b, MAX_SHIFT, dt)
        if np.isfinite(lag):
            lag += (j0 - (i0 - pad)) * dt      # correct for clipping at the array ends
        recs.append(dict(block=bi, t_mid=float(0.5 * (seg[0]["Start Time"] + seg[-1]["End Time"])),
                         n_lines=len(seg), lag_s=lag, peak=peak, margin=margin,
                         dur_med=float(np.median([d["End Time"] - d["Start Time"] for d in seg]))))
    if len(recs) < 4:
        return None, None
    df = pd.DataFrame(recs)
    good = df[(df.peak > PEAK_MIN) & (df.margin > MARGIN_MIN)].copy()
    if len(good) < 4:
        print(f"   [{sample}] only {len(good)}/{len(df)} blocks pass "
              f"(peak max {df.peak.max():.2f}, margin max {df.margin.max():.3f})")
        return None, df
    x = good.t_mid.values - ls[0]
    y = good.lag_s.values
    # Robust fit. A block whose correlation peak lands on the wrong lobe gives a wild
    # outlier (tens of seconds) while the rest sit on a line with tens of milliseconds
    # of scatter. A single rejection pass against an OLS fit is NOT safe here: one wild
    # point drags the fit far enough that good blocks land outside the band and get
    # discarded as collateral damage (found in review; it shifted one run by 6 ppm).
    # Remove the worst point one at a time instead, refitting after each removal, so
    # only genuine outliers ever leave.
    p = np.polyfit(x, y, 1)
    n_rej = 0
    while len(x) > 4:
        r = np.abs(y - np.polyval(p, x))
        w = int(np.argmax(r))
        if r[w] <= reject_tol:
            break
        x = np.delete(x, w)
        y = np.delete(y, w)
        n_rej += 1
        p = np.polyfit(x, y, 1)
    resid = y - np.polyval(p, x)
    # standard error of the slope, so runs can be compared against their own precision
    sxx = float(((x - x.mean()) ** 2).sum())
    dof = max(len(x) - 2, 1)
    slope_se = float(np.sqrt((resid**2).sum() / dof / sxx)) if sxx > 0 else np.nan
    speed = float(L[0].get("Scan Speed") or 0)
    spot = float(L[0].get("Spot Width") or 15)
    summary = dict(
        sample=sample.strip(),  # the io4 name can carry trailing whitespace; keep outputs clean
        n_lines=len(L),
        duration_h=round(float(ls[-1] - ls[0]) / 3600, 2),
        n_blocks_used=len(x),
        n_blocks_gate_failed=int(len(df) - len(good)),
        n_blocks_rejected_outlier=n_rej,
        slope_se_ppm=round(slope_se * 1e6, 1),
        lag_first_block_s=round(float(y[0]), 2),
        lag_last_block_s=round(float(y[-1]), 2),
        lag_change_s=round(float(y[-1] - y[0]), 2),
        drift_ppm=round(float(p[0]) * 1e6, 1),
        fit_resid_sd_s=round(float(resid.std()), 3),
        peak_min=round(float(good.peak.min()), 3),
        margin_min=round(float(good.margin.min()), 3),
        px_shift_end=round(float(p[0]) * float(ls[-1] - ls[0]) * speed / spot, 1),
    )
    df["sample"] = sample.strip()
    return summary, df
