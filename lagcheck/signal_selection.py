"""Optional choice of the matching signal by starting-alignment contrast.

For the CSV signal formats the default matching signal is chosen with
``ms.signal_selection = "highest_alignment_peak"``: the sum of the configured
columns and each configured column are compared by the peak of the automatic
starting-alignment search around the input placement, and the candidate with
the highest peak is used for the whole diagnosis. ``"sum"`` (the default for
netCDF input) uses the sum of the configured columns without comparison. The
search is the automatic one even when a starting alignment is declared, and a
candidate can win on peak height without meeting the starting-alignment
criteria; both facts are recorded. The choice is a declared one made from the
same record, not a verified pairing.
"""
import math

from .estimator import _coarse
from .models import InputError

MODES = ('sum', 'highest_alignment_peak')
TIE_TOLERANCE = 1e-9


def select_signal(summed, laser, analysis, load_columns):
    """Return (signal, record) for the candidate with the highest alignment peak.

    ``summed`` is the signal already loaded for all configured columns and
    ``load_columns(columns)`` loads the signal for a subset of them. The
    configured sum is kept unless a single column exceeds the best peak so far
    by more than ``TIE_TOLERANCE``, so equal candidates do not depend on
    floating-point rounding.
    """
    columns = list(summed.metadata['signal_columns'])
    window = analysis.get('coarse_window_s', 300.)
    if isinstance(window, bool) or not isinstance(window, (int, float)) or not math.isfinite(window) or not window > 0:
        raise InputError('coarse_window_s must be a positive finite number')
    candidates = [(columns, summed)]
    if len(columns) > 1:
        candidates += [([name], load_columns([name])) for name in columns]
    rows = []
    for cols, data in candidates:
        coarse = _coarse(data.time, data.signal, laser.lines, float(window))
        rows.append(dict(columns=cols, peak=coarse.get('peak'), margin=coarse.get('margin'),
                         offset_s=coarse.get('offset_s'), at_window_boundary=coarse.get('boundary'),
                         starting_alignment_criteria_met=bool(coarse.get('ok'))))
    best = 0
    for i, row in enumerate(rows):
        if row['peak'] is None:
            continue
        if rows[best]['peak'] is None or row['peak'] > rows[best]['peak'] + TIE_TOLERANCE:
            best = i
    record = dict(mode='highest_alignment_peak', configured_columns=columns,
                  selected_columns=rows[best]['columns'], candidates=rows, tie_tolerance=TIE_TOLERANCE,
                  criterion='Highest peak of the automatic starting-alignment search around the input placement',
                  note='Declared choice made from the same record; not an independent pairing. '
                       'The search does not use a declared starting alignment.')
    if all(row['peak'] is None for row in rows):
        record['note'] += ' No candidate showed temporal contrast; the configured signal was kept.'
    return candidates[best][1], record
