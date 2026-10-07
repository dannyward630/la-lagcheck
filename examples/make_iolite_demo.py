"""Create a public synthetic ESL example for the two-stage export (no iolite run)."""
import argparse
import csv
from datetime import datetime, timedelta
import json
from pathlib import Path

import numpy as np


def make_demo(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(148)
    intervals, cursor = [], 5.
    for _ in range(240):
        length = round(float(rng.uniform(5, 18)), 1)
        intervals.append((cursor, round(cursor + length, 1)))
        cursor = round(cursor + length + float(rng.uniform(3, 10)), 1)
    time = np.arange(0, cursor + 30, .1)
    signal = np.full(len(time), 2.)
    for start, end in intervals:
        signal[round((start + 3) * 10):round((end + 3) * 10)] = 1000.
    with (folder / 'signal.csv').open('x', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_s', 'signal'])
        writer.writerows(zip(time, signal))
    base = datetime(2026, 1, 1)

    def stamp(seconds):
        return (base + timedelta(seconds=seconds)).strftime('%Y-%m-%d %H:%M:%S.%f')

    with (folder / 'laser.csv').open('x', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['Timestamp', 'Comment', 'Sequence Number', 'Scan Velocity (um/s)', 'Laser State', 'X'])
        writer.writerow([stamp(0), 'Synthetic A', '1', '10', 'Off', 'metadata'])
        for index, (start, end) in enumerate(intervals):
            writer.writerow([stamp(start), '', '', '', 'On', str(index)])
            writer.writerow([stamp(end), '', '', '', 'Off', str(index)])
    config = dict(ms=dict(path='signal.csv', format='csv', header_row=0, encoding='utf-8',
                         delimiter=',', time_kind='numeric', time_unit='s', time_column='time_s',
                         signal_columns=['signal']),
                  laser=dict(path='laser.csv', format='laser_state_csv', header_row=0,
                             encoding='utf-8', delimiter=',', time_kind='datetime',
                             sample_selector='Synthetic A'),
                  analysis=dict(tolerance_s=.2, alignment_offset_s=3.))
    (folder / 'config.json').write_text(json.dumps(config, indent=2) + '\n', encoding='utf-8')
    (folder / 'truth.json').write_text(json.dumps(dict(
        rate_ppm=0., physical_offset_s=3., n_intervals=240, seed=148,
        ms_first_sample_at_laser_time=stamp(0),
        scope='Synthetic ESL formatting and fixed-coordinate timing control; not importer validation'),
        indent=2) + '\n', encoding='utf-8')
    return folder.resolve()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(make_demo(parser.parse_args().output))
