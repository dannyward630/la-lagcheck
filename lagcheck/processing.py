"""Transparent record segmentation; never edit, sort or join the source CSV."""
import csv
import hashlib
import io
from pathlib import Path

from .laser import _selected_rows, _column, _time
from .models import InputError
from .new_wave import FORMAT


def inspect_segments(path, config):
    path = Path(path)
    raw = path.read_bytes()
    text = raw.decode(config.get('encoding', 'utf-8-sig'))
    start = config.get('header_row', 0)
    if isinstance(start, bool) or not isinstance(start, int) or start < 0:
        raise InputError('header_row must be a nonnegative integer')
    try:
        reader = csv.reader(io.StringIO(''.join(text.splitlines(True)[start:])),
                            delimiter=config.get('delimiter', ','), strict=True)
        header = next(reader)
        rows = list(reader)
    except (StopIteration, csv.Error) as exc:
        raise InputError('Cannot inspect laser CSV: ' + str(exc)) from exc
    fmt = config.get('format', 'laser_state_csv')
    selected, labels, sample = _selected_rows(rows, header, dict(config, _inspect_segments=True), fmt)
    if not selected:
        raise InputError('No selected laser records')
    seq_col = next((i for i, h in enumerate(header) if h.strip() == 'Sequence Number'), None)
    sequences, inherited = [], ''
    for row in rows:
        if seq_col is not None:
            if seq_col >= len(row):
                raise InputError('Missing Sequence Number cell')
            inherited = row[seq_col].strip() or inherited
        sequences.append(inherited)
    groups = []
    for i in selected:
        if not groups or (fmt != FORMAT and (i != groups[-1][-1] + 1 or sequences[i] != sequences[groups[-1][-1]])):
            groups.append([])
        groups[-1].append(i)
    col = _column(header, config.get('time_column', 'Timestamp') if fmt in ('laser_state_csv', FORMAT) else config['start_column'])
    origin = _time(rows[selected[0]][col], config)[0]
    segments = []
    for number, group in enumerate(groups, 1):
        segment_origin = _time(rows[group[0]][col], config)[0]
        segment = dict(segment=number, start_record=group[0]+1, end_record=group[-1]+1,
                       sample=sample, sequence_id=sequences[group[0]] or None,
                       source_sha256=hashlib.sha256(raw).hexdigest(),
                       time_origin_seconds=str(segment_origin), origin_shift_s=float(segment_origin-origin),
                       source_record_count=len(group))
        if fmt in ('laser_state_csv', FORMAT):
            state_col = _column(header, config.get('state_column','Laser State'))
            segment['on_record_count'] = sum(rows[i][state_col].strip() == config.get('on_value','On') for i in group)
        if segments:
            segment['records_between_previous_segment'] = [segments[-1]['end_record']+1, group[0]]
            segment['gap_from_previous_last_record_s'] = float(segment_origin-_time(rows[groups[number-2][-1]][col], config)[0])
        segments.append(segment)
    return dict(version='0.4.0', source_path=str(path.resolve()), source_sha256=hashlib.sha256(raw).hexdigest(),
                original_selected_origin_seconds=str(origin), segments=segments,
                policy=('New Wave line/spot sequence changes remain inside the explicitly selected acquisition; explicit sample selection retains source gaps without bridging firing intervals.' if fmt == FORMAT else 'Each contiguous sample/sequence fragment is processed separately; no temporal bridge or combined correction.'))
