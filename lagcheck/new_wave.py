"""Explicit New Wave event-log semantics, without rewriting source records."""
from .models import InputError


FORMAT = 'new_wave_laser_log'


def sequence_labels(rows, header):
    columns = [i for i, h in enumerate(header) if h.strip() == 'Sequence Number']
    if len(columns) != 1:
        raise InputError('New Wave logs require one Sequence Number column')
    col = columns[0]
    labels, current = [], ''
    for row in rows:
        if col >= len(row):
            raise InputError('Missing New Wave Sequence Number cell')
        current = row[col].strip() or current
        labels.append(current)
    return labels


def identity(rows, header, indices, sample, config):
    if config.get('expected_sequence') is not None:
        raise InputError('New Wave Sequence Number identifies a line/spot, not an acquisition; do not set expected_sequence')
    sequences = sequence_labels(rows, header)
    selected = [sequences[i] for i in indices]
    if any(not v for v in selected):
        raise InputError('New Wave selected records have an undefined line/spot Sequence Number')
    return dict(sample=sample, sequence_id=None, sequence_column='Sequence Number',
                sequence_verified=False, sequence_ids=list(dict.fromkeys(selected)),
                scope='One explicitly paired input file; Sequence Number identifies scan lines or spots, not acquisitions',
                adapter='new_wave_events_v1', sample_scope=config.get('sample_scope', 'selected_sample'))


def annotate(laser, source_indices, labels, sequences):
    table = laser.table
    header = [h.strip() for h in table['header']]
    fields = ['X(um)', 'Y(um)', 'Vertix Number', 'Scan Velocity (um/s)']
    columns = {name: header.index(name) for name in fields if name in header}
    records = []
    updates = 0
    for number, (line, (a, b)) in enumerate(zip(laser.lines, table['interval_rows']), 1):
        line['Sample Name'] = labels[source_indices[a]]
        vertices = []
        for i in range(a, b + 1):
            row = table['rows'][i]
            vertices.append(dict(source_record=source_indices[i] + 1,
                time_s=float(table['times'][i][0] - table['origin']),
                fields={name: row[col] if col < len(row) else '' for name, col in columns.items()}))
        updates += max(0, b - a - 1)
        records.append(dict(interval=number, sample=line['Sample Name'],
            sequence_id=sequences[source_indices[a]], start_record=source_indices[a] + 1,
            end_record=source_indices[b] + 1, start_s=line['Start Time'], end_s=line['End Time'],
            event_records=vertices))
    laser.metadata['new_wave'] = dict(adapter='new_wave_events_v1', intervals=records,
        interval_count=len(records), state_update_count=updates,
        rule='First On to first subsequent Off within the same sample and Sequence Number; intermediate On records are state/vertex updates, not new firing intervals.',
        geometry_note='Recorded X/Y and vertices are retained; reconstruction assumptions are stated in the view metadata.')
    laser.metadata['input_processing']['selection_record_gaps'] = [dict(after_record=a+1, before_record=b+1, omitted_record_count=b-a-1) for a,b in zip(source_indices,source_indices[1:]) if b != a+1]
    laser.metadata['input_processing']['event_semantics'] = laser.metadata['new_wave']['rule']
    return laser
