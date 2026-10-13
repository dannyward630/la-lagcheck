"""Tolerance requests in seconds or micrometres, and distances equivalent to time differences.

A request in micrometres is converted to seconds once, before the diagnosis.
The estimator receives seconds only, so the same seconds give the same result
however they were requested.
"""
import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from lagcheck.cli import run as _run, main, digest
from lagcheck.estimator import diagnose
from lagcheck.iolite import export_iolite
from lagcheck.laser import load_laser
from lagcheck.models import InputError, LaserData
from lagcheck.preview_geometry import DEFAULT_PITCH_UM
from lagcheck.tolerance import (estimator_settings, firing_speed, read_request, resolve, spot_width, distance_scale,
                                DISTANCE_RULE, SECONDS_RULE)
from examples.make_iolite_demo import make_demo
import test_cli
from test_estimator import synthetic

HEADER = ['Timestamp', 'Comment', 'Sequence Number', 'Scan Velocity (um/s)', 'Laser State', 'X(um)', 'Y(um)',
          'Spot Type', 'Spot Size (um)', 'Spot Angle (deg)']
STATE = dict(format='laser_state_csv', time_kind='numeric', time_unit='s', sample_selector='Synthetic A')
SPEED_SOURCE = "laser log column 'Scan Velocity (um/s)' at each firing start; a blank cell takes the last value above it"


def run(*args, **kwargs):
    kwargs.setdefault('preview', False)
    return _run(*args, **kwargs)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def state_log(path, intervals, speeds=375, positioning=20000, spot=('', '15 x 15', ''), move=(1., 0.), header=HEADER, on_speeds=None):
    """An ESL-style state log with numeric times.

    Before each firing interval two Off records position the stage: the first
    carries the positioning speed, the second the scan speed. The On and Off
    records of the firing leave the speed cell blank, unless ``on_speeds``
    writes a speed in the On record itself. ``speeds`` and ``on_speeds`` are
    one value, or one value for each interval; a ``speeds`` of None leaves out
    the record that sets it. ``move`` is the direction of the scan, (0, 0) for
    a spot.
    """
    speeds = list(speeds) if isinstance(speeds, (list, tuple)) else [speeds]*len(intervals)
    on_speeds = list(on_speeds) if isinstance(on_speeds, (list, tuple)) else [on_speeds]*len(intervals)
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(header)
        cells = dict(zip(HEADER[7:], spot))
        w.writerow([c for name, c in zip(HEADER, ['0.000', 'Synthetic A', '1', '', 'Off', '0.0', '0.0', *spot]) if name in header])
        moving = 100.
        for i, ((a, b), speed, on_speed) in enumerate(zip(intervals, speeds, on_speeds)):
            dx, dy = move[i] if isinstance(move, list) else move
            for value in (speed, on_speed):
                if isinstance(value, (int, float)) and value:
                    moving = abs(float(value))                    # the stage keeps the last speed that was set
            length = moving*(b-a)
            forward = i % 2 == 0                                  # serpentine: every second line runs backwards
            start = (0., 20.*i) if forward else (dx*length, 20.*i+dy*length)
            end = (dx*length, 20.*i+dy*length) if forward else (0., 20.*i)
            rows = [[f'{a-.02:.3f}', '', '', positioning, 'Off', *start], [f'{a-.01:.3f}', '', '', speed, 'Off', *start],
                    [f'{a:.3f}', '', '', '' if on_speed is None else on_speed, 'On', *start], [f'{b:.3f}', '', '', '', 'Off', *end]]
            if positioning is None:
                del rows[0]
            if speed is None:
                rows = [row for row in rows if row[3] != speed or row[4] != 'Off' or row[0] == f'{b:.3f}']
            for row in rows:
                full = dict(zip(HEADER[:7], [row[0], row[1], row[2], '' if row[3] is None else row[3], row[4], f'{row[5]:.1f}', f'{row[6]:.1f}']), **cells)
                w.writerow([full[name] for name in header])
    return Path(path)


def small_log(path, speeds=375, n=3, **kwargs):
    return load_laser(state_log(path, [(10.+20*i, 20.+20*i) for i in range(n)], speeds, **kwargs), STATE)


def rewrite(path, change):
    """Pass the records of a written log, header included, through ``change``, which edits the list in place."""
    with Path(path).open(encoding='utf-8', newline='') as f:
        rows = list(csv.reader(f))
    change(rows)
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        csv.writer(f).writerows(rows)
    return Path(path)


def blank_first_on_x(rows):
    """Remove one coordinate: the x position of the first On record."""
    next(row for row in rows[1:] if row[4] == 'On')[5] = ''


def write_scan_case(root, speeds=375, spot=('', '15 x 15', ''), move=(1., 0.), **analysis):
    """Signal and state log of 240 scan lines with a rate of 100 µs/s, and a configuration for them."""
    root.mkdir(parents=True, exist_ok=True)
    ms, laser = synthetic(100., 3.)
    with (root/'signal.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_ms', 'element'])
        w.writerows(zip(ms.time*1000, ms.signal))
    state_log(root/'laser.csv', [(row['Start Time'], row['End Time']) for row in laser.lines], speeds, spot=spot, move=move)
    config = {'ms': {'path': 'signal.csv', 'format': 'csv', 'header_row': 0, 'encoding': 'utf-8', 'delimiter': ',',
                     'time_kind': 'numeric', 'time_unit': 'ms', 'time_column': 'time_ms', 'signal_columns': ['element']},
              'laser': dict(STATE, path='laser.csv'),
              'analysis': dict({'alignment_offset_s': 3., 'skip_record_edge_lines': False}, **analysis)}
    (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
    return config


class Case(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def scan(self, name, **kwargs):
        write_scan_case(self.root/name, **kwargs)
        return self.root/name

    def interval(self, name, **analysis):
        root = self.root/name
        root.mkdir()
        config = test_cli.WorkflowTests.make_case(None, root)
        config['analysis'] = dict({'alignment_offset_s': 8., 'skip_record_edge_lines': False}, **analysis)
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        return root

    def command(self, *args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = main([str(a) for a in args])
        return code, out.getvalue(), err.getvalue()

    def error(self, *args):
        code, printed, errors = self.command(*args)
        self.assertEqual((code, printed), (1, ''))
        return json.loads(errors.splitlines()[-1])['error']


class SpeedTests(Case):
    def test_one_speed_converts_micrometres_to_the_seconds_of_the_study_examples(self):
        for distance, speed in ((75., 375), (45., 225)):
            laser = small_log(self.root/f'{speed}.csv', speed)
            settings, record, evidence = resolve(read_request({'tolerance_um': distance}), {'tolerance_um': distance}, laser)
            self.assertEqual(settings, {'tolerance_s': .2})
            self.assertEqual(record, dict(requested_unit='um', requested_value=distance, request_source='analysis.tolerance_um',
                                          overridden_configuration_value=None, resolved_tolerance_s=.2, rule=DISTANCE_RULE,
                                          speed_um_s=float(speed), speed_source=SPEED_SOURCE, reused_from_saved_configuration=False))
            self.assertEqual((evidence['status'], evidence['min_um_s'], evidence['max_um_s'], evidence['values_um_s'], evidence['n_firing_intervals']),
                             ('credible', float(speed), float(speed), [float(speed)], 3))

    def test_several_speeds_use_the_largest_and_report_the_range(self):
        laser = small_log(self.root/'mixed.csv', [225, 375, 225])
        evidence = firing_speed(laser)
        self.assertEqual((evidence['status'], evidence['min_um_s'], evidence['max_um_s'], evidence['values_um_s']), ('credible', 225., 375., [225., 375.]))
        settings, record, _ = resolve(read_request({'tolerance_um': 75.}), {'tolerance_um': 75.}, laser)
        self.assertEqual((settings['tolerance_s'], record['speed_um_s']), (.2, 375.))       # not 75/225 = 0.33 s
        scale = distance_scale({}, [], .2, 100., evidence, spot_width(laser, {}))
        row = scale['quantities'][0]
        self.assertEqual((row['quantity'], row['seconds']), ('tolerance', .2))
        self.assertAlmostEqual(row['um_min'], 45.)
        self.assertAlmostEqual(row['um_max'], 75.)

    def test_positioning_speed_of_an_off_record_is_not_a_firing_speed(self):
        laser = small_log(self.root/'log.csv', 375, positioning=20000)
        column = [row[3] for row in laser.table['rows']]
        self.assertIn('20000', column)                               # the log holds the positioning speed
        evidence = firing_speed(laser)
        self.assertEqual((evidence['values_um_s'], evidence['max_um_s'], evidence['log_values_um_s']), ([375.], 375., [375.]))
        # Without the record that sets the scan speed, the firing start takes the positioning speed from the record above it.
        # The recorded positions and times contradict that speed, so the interval has no firing speed.
        inherited = small_log(self.root/'inherited.csv', [375, None, 375])
        self.assertEqual([row[3] for row in inherited.table['rows'] if row[4] == 'Off' and row[3]], ['20000', '375', '20000', '20000', '375'])
        evidence = firing_speed(inherited)
        self.assertEqual((evidence['status'], evidence['log_values_um_s'], evidence['reason'][:6]), ('not_computable', [375.], '1 of 3'))
        # Without positions nothing in the log tells that speed from a scan speed, so no speed taken from a record above is used.
        without_positions = small_log(self.root/'no_positions.csv', [375, None, 375], header=HEADER[:5]+HEADER[7:])
        evidence = firing_speed(without_positions)
        self.assertEqual((evidence['status'], evidence['values_um_s'], evidence['log_values_um_s'], evidence['reason'][:6]),
                         ('not_computable', [], [], '3 of 3'))
        # A speed written in the Off record that ends a firing belongs to the next move.
        path = state_log(self.root/'after.csv', [(10., 20.), (30., 40.)], 375)
        rows = list(csv.reader(path.open(encoding='utf-8')))
        for row in rows:
            if row[4] == 'Off' and row[0] in ('20.000', '40.000'):
                row[3] = '20000'
        with path.open('w', newline='', encoding='utf-8') as f:
            csv.writer(f).writerows(rows)
        self.assertEqual(firing_speed(load_laser(path, STATE))['values_um_s'], [375.])

    def test_packaged_logs_give_the_scan_speed_and_not_the_positioning_speed(self):
        package = Path(__file__).resolve().parents[1]
        grating = read(package/'examples/chromium_grating/auto.json')['laser']
        logs = [(package/'examples/R01/R01_native_laser_log.csv', dict(format='laser_state_csv', time_kind='numeric', time_unit='s'), 903),
                (package/'examples/chromium_grating'/grating['path'], {k: v for k, v in grating.items() if k != 'path'}, 621)]
        for path, config, intervals in logs:
            laser = load_laser(path, dict(config, zero_duration_policy='exclude'))
            header = [name.strip() for name in laser.table['header']]
            recorded = {row[header.index('Scan Velocity (um/s)')].strip() for row in laser.table['rows']}
            self.assertEqual(recorded, {'', '375', '20000'}, path.name)
            evidence = firing_speed(laser)
            self.assertEqual((evidence['status'], evidence['values_um_s'], evidence['n_firing_intervals']), ('credible', [375.], intervals), path.name)
            settings, record, _ = resolve(read_request({'tolerance_um': 75.}), {'tolerance_um': 75.}, laser)
            self.assertEqual((settings['tolerance_s'], record['speed_um_s']), (.2, 375.), path.name)
            self.assertEqual(spot_width(laser, {})['status'], 'not_computable', path.name)       # '15 x 15' without a recorded spot type

    def test_direction_of_the_scan_does_not_change_the_speed(self):
        backward = small_log(self.root/'signed.csv', [375, -375, 375])        # a log that signs the speed of a return line
        self.assertEqual(firing_speed(backward)['values_um_s'], [375.])
        for name, move in (('x', (1., 0.)), ('y', (0., 1.)), ('diagonal', (.6, .8))):
            evidence = firing_speed(small_log(self.root/f'{name}.csv', 375, move=move))
            self.assertEqual((evidence['status'], evidence['max_um_s']), ('credible', 375.), name)
        row = distance_scale({}, [], .2, 100., evidence, dict(status='not_computable', width_um=None, source=None, reason='x'))['quantities'][0]
        self.assertEqual(sorted(row), ['basis', 'limit_s', 'quantity', 'seconds', 'spot_widths_max', 'spot_widths_min', 'um_max', 'um_min'])
        self.assertAlmostEqual(row['um_max'], 75.)                    # a length along the scan, whatever its direction

    def test_unresolved_speed_of_an_interval_is_not_replaced_by_the_first_one(self):
        laser = small_log(self.root/'gap.csv', [375, 'fast', 375])       # the second line has no readable speed
        evidence = firing_speed(laser)
        self.assertEqual((evidence['status'], evidence['max_um_s'], evidence['log_values_um_s']), ('not_computable', None, [375.]))
        self.assertEqual(evidence['reason'], "1 of 3 firing intervals have no usable scan speed in column 'Scan Velocity (um/s)' "
                         '(1 magnitude not a positive finite number) and analysis.scan_speed_um_s is not set')
        with self.assertRaisesRegex(InputError, 'cannot be converted to seconds for this record: 1 of 3 firing intervals'):
            resolve(read_request({'tolerance_um': 75.}), {'tolerance_um': 75.}, laser)
        for bad in (0, 'inf', 'nan'):
            self.assertEqual(firing_speed(small_log(self.root/f'bad_{bad}.csv', [375, bad, 375], positioning=None))['status'], 'not_computable', bad)
        # A blank cell is not a missing speed: it takes the last value above it. Without a value above, the interval has none.
        self.assertEqual(firing_speed(small_log(self.root/'blank.csv', [375, '', 375], positioning=None))['values_um_s'], [375.])
        self.assertEqual(firing_speed(small_log(self.root/'first.csv', ['', 375, 375], positioning=None))['reason'][:6], '1 of 3')
        # A speed that changes while the laser fires leaves the speed of that interval undefined.
        path = state_log(self.root/'change.csv', [(10., 20.), (30., 40.)], 375)
        rows = list(csv.reader(path.open(encoding='utf-8')))
        on = next(i for i, row in enumerate(rows) if row[4] == 'On')
        rows.insert(on+1, ['15.000', '', '', '100', '', '0.0', '0.0', '', '15 x 15', ''])       # a record without a state, during the firing
        with path.open('w', newline='', encoding='utf-8') as f:
            csv.writer(f).writerows(rows)
        self.assertEqual(firing_speed(load_laser(path, STATE))['reason'][:6], '1 of 2')

    def test_positions_are_kept_for_each_interval_and_a_missing_one_hides_no_other(self):
        unused = dict.fromkeys(('not_recorded', 'not_a_positive_finite_number', 'changed_during_firing', 'duration_not_positive',
                                'contradicted_by_positions', 'inherited_without_positions'), 0)
        inherited = dict(written_at_firing_start=0, inherited_from_record_above=3, not_recorded=0)
        lines, request = [(10., 20.), (30., 40.), (50., 60.)], {'tolerance_um': 75.}
        # In each log the second firing start has no record that sets the scan speed, so it takes the positioning speed of 20000 µm/s.
        # 1. Every coordinate is recorded: the positions of the second interval contradict that speed.
        complete = small_log(self.root/'complete.csv', [375, None, 375])
        evidence = firing_speed(complete)
        self.assertEqual((evidence['status'], evidence['intervals_without_positions'], evidence['log_speed_origins']), ('not_computable', 0, inherited))
        self.assertEqual(evidence['log_speeds_not_used'], dict(unused, contradicted_by_positions=1))
        self.assertEqual(evidence['reason'], "1 of 3 firing intervals have no usable scan speed in column 'Scan Velocity (um/s)' "
                         '(1 contradicted by the recorded positions and times by more than a factor of two) and analysis.scan_speed_um_s is not set')
        # 2. The log has no position columns: no interval has a position, and no speed taken from a record above is used.
        absent = small_log(self.root/'absent.csv', [375, None, 375], header=HEADER[:5]+HEADER[7:])
        evidence = firing_speed(absent)
        self.assertEqual((evidence['status'], evidence['intervals_without_positions'], evidence['max_um_s']), ('not_computable', 3, None))
        self.assertEqual(evidence['log_speeds_not_used'], dict(unused, inherited_without_positions=3))
        self.assertEqual(evidence['reason'], "3 of 3 firing intervals have no usable scan speed in column 'Scan Velocity (um/s)' "
                         '(3 taken from a record above the firing start, with no position of the interval to check it) '
                         'and analysis.scan_speed_um_s is not set')
        # 3. One coordinate of the first interval is missing: the second interval is still checked against its own positions.
        partial = load_laser(rewrite(state_log(self.root/'partial.csv', lines, [375, None, 375]), blank_first_on_x), STATE)
        evidence = firing_speed(partial)
        self.assertEqual((evidence['status'], evidence['intervals_without_positions'], evidence['log_values_um_s']), ('not_computable', 1, [375.]))
        self.assertEqual(evidence['log_speeds_not_used'], dict(unused, contradicted_by_positions=1, inherited_without_positions=1))
        self.assertEqual(evidence['reason'], "2 of 3 firing intervals have no usable scan speed in column 'Scan Velocity (um/s)' "
                         '(1 contradicted by the recorded positions and times by more than a factor of two; '
                         '1 taken from a record above the firing start, with no position of the interval to check it) '
                         'and analysis.scan_speed_um_s is not set')
        for laser in (complete, absent, partial):
            # 75 µm is not converted with the positioning speed (75/20000 = 0.00375 s): the request is an input error.
            with self.assertRaises(InputError) as caught:
                resolve(read_request(request), request, laser)
            self.assertTrue(str(caught.exception).endswith('Give the tolerance in seconds (analysis.tolerance_s or --tolerance-s), '
                                                           'or declare analysis.scan_speed_um_s.'), str(caught.exception))
            # With the scan speed declared it is 75/375 = 0.2 s, and a tolerance in seconds needs no speed at all.
            declared = dict(request, scan_speed_um_s=375.)
            settings, record, evidence = resolve(read_request(declared), declared, laser)
            self.assertEqual((settings['tolerance_s'], record['speed_source'], evidence['status']), (.2, 'analysis.scan_speed_um_s (declared)', 'credible'))
            self.assertEqual(resolve(read_request({'tolerance_s': .2}), {'tolerance_s': .2}, laser)[0], {'tolerance_s': .2})
        # 4. A missing coordinate in the first interval and a spot in the second: the record is a spot record, as it is without the gap.
        path = state_log(self.root/'spot.csv', lines, 375, move=[(1., 0.), (0., 0.), (1., 0.)])
        spot = load_laser(rewrite(path, blank_first_on_x), STATE)
        for speed in (None, 375.):
            evidence = firing_speed(spot, speed)
            self.assertEqual((evidence['status'], evidence['declared_speed_accepted'], evidence['intervals_without_positions']), ('not_computable', False, 1))
            self.assertTrue(evidence['reason'].startswith('1 of 3 firing intervals show no stage movement (spot record)'), evidence['reason'])
        declared = dict(request, scan_speed_um_s=375.)
        with self.assertRaises(InputError) as caught:
            resolve(read_request(declared), declared, spot)
        self.assertTrue(str(caught.exception).endswith('Give the tolerance in seconds (analysis.tolerance_s or --tolerance-s).'), str(caught.exception))
        self.assertEqual(resolve(read_request({'tolerance_s': .2}), {'tolerance_s': .2}, spot)[0], {'tolerance_s': .2})

    def test_speed_written_in_the_on_record_is_used_and_one_taken_from_above_needs_the_positions_of_its_interval(self):
        lines, request = [(10., 20.), (30., 40.), (50., 60.)], {'tolerance_um': 75.}
        written = dict(written_at_firing_start=3, inherited_from_record_above=0, not_recorded=0)
        # The On record itself states 375 µm/s, below a record with the positioning speed. It is used with positions,
        # without position columns, and when one interval has lost a coordinate.
        logs = [('positions', small_log(self.root/'on.csv', None, on_speeds=375), 0),
                ('no_columns', small_log(self.root/'on_plain.csv', None, on_speeds=375, header=HEADER[:5]+HEADER[7:]), 3),
                ('one_gap', load_laser(rewrite(state_log(self.root/'on_gap.csv', lines, None, on_speeds=375), blank_first_on_x), STATE), 1)]
        for name, laser, gaps in logs:
            self.assertIn('20000', [row[3] for row in laser.table['rows']], name)
            evidence = firing_speed(laser)
            self.assertEqual((evidence['status'], evidence['values_um_s'], evidence['source']), ('credible', [375.], SPEED_SOURCE), name)
            self.assertEqual((evidence['log_speed_origins'], evidence['intervals_without_positions']), (written, gaps), name)
            self.assertEqual(sum(evidence['log_speeds_not_used'].values()), 0, name)
            self.assertEqual(resolve(read_request(request), request, laser)[0]['tolerance_s'], .2, name)
        # Positions still set aside a speed that they contradict, wherever it was written.
        def positioning_speed_in_the_on_records(rows):
            for row in rows[1:]:
                if row[4] == 'On':
                    row[3] = '20000'
        wrong = load_laser(rewrite(state_log(self.root/'on_wrong.csv', lines, None, on_speeds=375), positioning_speed_in_the_on_records), STATE)
        evidence = firing_speed(wrong)
        self.assertEqual((evidence['status'], evidence['log_speeds_not_used']['contradicted_by_positions'], evidence['log_speed_origins']),
                         ('not_computable', 3, written))
        # A speed taken from the record above is used when the positions of its own interval agree with it.
        supported = small_log(self.root/'above.csv', 375)
        evidence = firing_speed(supported)
        self.assertEqual((evidence['status'], evidence['values_um_s'], evidence['log_speed_origins']['inherited_from_record_above']), ('credible', [375.], 3))
        self.assertEqual(resolve(read_request(request), request, supported)[0]['tolerance_s'], .2)
        # Without the positions of one interval, that interval has no speed, although the two others agree with theirs.
        gap = load_laser(rewrite(state_log(self.root/'above_gap.csv', lines, 375), blank_first_on_x), STATE)
        evidence = firing_speed(gap)
        self.assertEqual((evidence['status'], evidence['log_values_um_s'], evidence['log_speeds_not_used']['inherited_without_positions']),
                         ('not_computable', [375.], 1))
        with self.assertRaisesRegex(InputError, '1 of 3 firing intervals have no usable scan speed'):
            resolve(read_request(request), request, gap)

    def test_scan_fields_that_cannot_be_read_do_not_stop_a_run_in_seconds(self):
        laser = small_log(self.root/'log.csv', 375)
        # Firing intervals that do not correspond to their records one for one cannot be given a speed or a position.
        shorter = LaserData(lines=laser.lines[:-1], metadata=laser.metadata, table=laser.table)
        evidence = firing_speed(shorter)
        self.assertEqual((evidence['status'], evidence['reason'], evidence['log_speed_origins'], evidence['intervals_without_positions']),
                         ('not_computable', 'the scan fields of the laser log could not be read', None, None))
        self.assertEqual(resolve(read_request({'tolerance_s': .2}), {'tolerance_s': .2}, shorter)[0], {'tolerance_s': .2})
        # The message asks for seconds only: a declared speed is not applied to a record whose scan fields cannot be read.
        with self.assertRaisesRegex(InputError, r'the scan fields of the laser log could not be read\. Give the tolerance in seconds \(analysis\.tolerance_s or --tolerance-s\)\.$'):
            resolve(read_request({'tolerance_um': 75.}), {'tolerance_um': 75.}, shorter)
        self.assertEqual((firing_speed(shorter, 375.)['status'], firing_speed(shorter, 375.)['declared_speed_accepted']), ('not_computable', False))
        # A firing interval without a positive duration has no speed.
        lines = [dict(line) for line in laser.lines]
        lines[1]['End Time'] = lines[1]['Start Time']
        evidence = firing_speed(LaserData(lines=lines, metadata=laser.metadata, table=laser.table))
        self.assertEqual((evidence['status'], evidence['log_speeds_not_used']['duration_not_positive'], evidence['reason'][:6]), ('not_computable', 1, '1 of 3'))

    def test_declared_speed_is_used_and_a_difference_from_the_log_is_recorded(self):
        laser = small_log(self.root/'log.csv', 375)
        evidence = firing_speed(laser, 300.)
        self.assertEqual((evidence['status'], evidence['source'], evidence['values_um_s'], evidence['log_values_um_s'], evidence['declared_differs_from_log']),
                         ('credible', 'analysis.scan_speed_um_s (declared)', [300.], [375.], True))
        analysis = {'tolerance_um': 75., 'scan_speed_um_s': 300.}
        settings, record, _ = resolve(read_request(analysis), analysis, laser)
        self.assertEqual((settings, record['speed_um_s'], record['speed_source']),
                         ({'scan_speed_um_s': 300., 'tolerance_s': .25}, 300., 'analysis.scan_speed_um_s (declared)'))
        self.assertIs(firing_speed(laser, 375.)['declared_differs_from_log'], False)
        without_column = small_log(self.root/'plain.csv', 375, header=HEADER[:3]+HEADER[4:])
        self.assertEqual(firing_speed(without_column)['reason'], 'the laser log gives no scan speed and analysis.scan_speed_um_s is not set')
        self.assertEqual(firing_speed(without_column, 375.)['declared_differs_from_log'], None)
        for invalid in (0., -5., True, float('nan'), '375'):          # not a usable declaration; the estimator reports it as before
            self.assertEqual(firing_speed(without_column, invalid)['status'], 'not_computable', invalid)

    def test_spot_records_and_new_wave_logs_have_no_distance_and_need_seconds(self):
        spots = small_log(self.root/'spots.csv', 375, move=(0., 0.))
        mixed = small_log(self.root/'mixed.csv', 375, move=[(1., 0.), (0., 0.), (1., 0.)])
        for laser, count in ((spots, 3), (mixed, 1)):
            for declared in (None, 375.):
                evidence = firing_speed(laser, declared)
                self.assertEqual((evidence['status'], evidence['max_um_s'], evidence['declared_speed_accepted']), ('not_computable', None, False))
                self.assertTrue(evidence['reason'].startswith(f'{count} of 3 firing intervals show no stage movement (spot record)'))
            with self.assertRaises(InputError) as caught:
                resolve(read_request({'tolerance_um': 75.}), {'tolerance_um': 75., 'scan_speed_um_s': 375.}, laser)
            self.assertTrue(str(caught.exception).endswith('Give the tolerance in seconds (analysis.tolerance_s or --tolerance-s).'))
            settings, record, _ = resolve(read_request({'tolerance_s': .2}), {'tolerance_s': .2}, laser)
            self.assertEqual((settings, record['rule'], record['resolved_tolerance_s']), ({'tolerance_s': .2}, SECONDS_RULE, .2))
        rows = [[0, '1', '1', '', 'Synthetic A', 0, 0, '', 'Off'], [1, '', '', '', '', 0, 0, '', 'On'],
                [1.1, '', '', '2', '', 0, 0, '20000', 'On'], [5, '', '', '', '', 4, 0, '', 'Off']]
        path = self.root/'new_wave.csv'
        with path.open('w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['Timestamp', 'Sequence Number', 'SubPoint Number', 'Vertix Number', 'Comment', 'X(um)', 'Y(um)',
                        'Scan Velocity (um/s)', 'Laser State'])
            w.writerows(rows)
        new_wave = load_laser(path, dict(format='new_wave_laser_log', time_kind='numeric', time_unit='s', sample_scope='all'))
        for declared in (None, 375.):
            evidence = firing_speed(new_wave, declared)
            self.assertEqual((evidence['status'], evidence['declared_speed_accepted']), ('not_computable', False))
            self.assertIn('no firing speed is taken from this format', evidence['reason'])
        with self.assertRaisesRegex(InputError, 'New Wave logs can carry a positioning speed'):
            resolve(read_request({'tolerance_um': 75.}), {'tolerance_um': 75., 'scan_speed_um_s': 375.}, new_wave)
        self.assertEqual(resolve(read_request({'tolerance_s': .2}), {'tolerance_s': .2}, new_wave)[0], {'tolerance_s': .2})


class SpotWidthTests(Case):
    def width(self, name, spot, analysis=None, **kwargs):
        return spot_width(small_log(self.root/f'{name}.csv', 375, spot=spot, **kwargs), analysis or {})

    def test_width_is_taken_from_the_log_only_when_shape_size_angle_and_direction_leave_no_choice(self):
        known = [('circle', ('Circle', '15', '')), ('circle_pair', ('circle', '15 x 15', '')), ('square', ('Square', '15 x 15', '0')),
                 ('square_quarter_turn', ('Rectangle', '15 x 15', '90'))]
        for name, spot in known:
            record = self.width(name, spot)
            self.assertEqual((record['status'], record['width_um'], record['reason']), ('known', 15., None), name)
            self.assertIn("laser log columns 'Spot Size (um)', 'Spot Type' and 'Spot Angle (deg)'", record['source'])
        unknown = [('blank_type', ('', '15 x 15', ''), 'the spot shape is not recorded'), ('turned', ('Square', '15 x 15', '45'), 'turned by 45 degrees'),
                   ('no_angle', ('Square', '15 x 15', ''), 'the spot angle is not recorded'), ('sides', ('Rectangle', '10 x 20', '0'), 'which side'),
                   ('ellipse', ('Ellipse', '10 x 20', '0'), "spot type 'Ellipse' is not interpreted"), ('text', ('Circle', 'large', ''), 'not read as'),
                   ('oval', ('Circle', '10 x 20', ''), 'two different sizes'), ('no_size', ('Circle', '', ''), 'not read as')]
        for name, spot, reason in unknown:
            record = self.width(name, spot)
            self.assertEqual((record['status'], record['width_um'], record['source']), ('not_computable', None, None), name)
            self.assertIn(reason, record['reason'], name)
        self.assertIn('not recorded as a stage axis', self.width('diagonal', ('Square', '15 x 15', '0'), move=(.6, .8))['reason'])
        self.assertEqual(self.width('y_scan', ('Square', '15 x 15', '0'), move=(0., 1.))['status'], 'known')

    def test_several_spot_sizes_or_a_log_without_the_field_give_no_width(self):
        path = state_log(self.root/'two.csv', [(10., 20.), (30., 40.)], 375, spot=('Circle', '15', ''))
        rows = list(csv.reader(path.open(encoding='utf-8')))
        for row in rows[len(rows)//2:]:
            row[8] = '30'
        with path.open('w', newline='', encoding='utf-8') as f:
            csv.writer(f).writerows(rows)
        self.assertEqual(spot_width(load_laser(path, STATE), {})['reason'], 'the firing intervals record more than one spot size, type or angle')
        without = small_log(self.root/'without.csv', 375, header=HEADER[:7])
        self.assertIn("no column 'Spot Size (um)' is read from this laser log", spot_width(without, {})['reason'])

    def test_size_or_angle_that_is_not_a_finite_number_gives_no_width(self):
        for angle in ('nan', 'NaN', 'inf', '-inf', 'Infinity', '1e999'):
            record = self.width(f'angle_{angle}', ('Square', '15 x 15', angle))
            self.assertEqual((record['status'], record['width_um'], record['source']), ('not_computable', None, None), angle)
            self.assertEqual(record['reason'], f'the spot angle {angle!r} is not a finite number', angle)
        digits = '9'*400                                             # a run of digits that overflows when it is read as a number
        for name, spot in (('circle', ('Circle', digits, '')), ('square', ('Square', f'{digits} x {digits}', '0'))):
            laser = small_log(self.root/f'overflow_{name}.csv', 375, spot=spot)
            record = spot_width(laser, {})
            self.assertEqual((record['status'], record['width_um'], record['source']), ('not_computable', None, None), name)
            self.assertTrue(record['reason'].endswith('more characters is not a finite number'), record['reason'])
            self.assertLess(len(record['reason']), 120, name)             # the message does not repeat the whole cell
            # The distance in micrometres is still given. A number of spot widths is not.
            row = distance_scale({}, [], .2, 100., firing_speed(laser), record)['quantities'][0]
            self.assertEqual((row['um_max'], row['spot_widths_min'], row['spot_widths_max']), (75., None, None), name)
        for size in ('0', '15 x 0', '0 x 15'):
            self.assertIn('is not read as one or two positive numbers', self.width(f'size_{size}', ('Circle', size, ''))['reason'], size)
        # Finite angles that leave the sides along the stage axes give the side as before, and so does a circle.
        for angle in ('0', '90', '180', '270', '-90', '360.0'):
            record = self.width(f'turn_{angle}', ('Square', '15 x 15', angle))
            self.assertEqual((record['status'], record['width_um'], record['reason']), ('known', 15., None), angle)
        record = self.width('circle', ('Circle', '15', ''))
        self.assertEqual((record['status'], record['width_um'], record['reason']), ('known', 15., None))
        # A width that is not a positive finite number never yields a number of spot widths, whatever supplied it.
        speed = firing_speed(small_log(self.root/'speed.csv', 375))
        for width in (float('inf'), float('nan'), 0., -15.):
            row = distance_scale({}, [], .2, 100., speed, dict(status='known', width_um=width, source='x', reason=None))['quantities'][0]
            self.assertEqual((row['um_max'], row['spot_widths_min'], row['spot_widths_max']), (75., None, None), width)

    def test_spot_that_changes_while_the_laser_fires_gives_no_width(self):
        def during_first_firing(cells):
            """Insert a record without a state after the first On record; ``cells`` are its spot type, size and angle."""
            def change(rows):
                on = next(i for i, row in enumerate(rows) if row[4] == 'On')
                rows.insert(on+1, ['15.000', '', '', '', '', rows[on][5], rows[on][6], *cells])
            return change
        lines = [(10., 20.), (30., 40.), (50., 60.)]
        # The size goes from 15 to 30 during the first firing and is 15 again in the record that ends it. Type and angle likewise.
        cases = [('size', ('Circle', '15', ''), ('', '30', '')), ('type', ('Circle', '15', ''), ('Square', '', '')),
                 ('angle', ('Square', '15 x 15', '0'), ('', '', '45'))]
        for name, spot, cells in cases:
            laser = load_laser(rewrite(state_log(self.root/f'{name}.csv', lines, 375, spot=spot), during_first_firing(cells)), STATE)
            header = [h.strip() for h in laser.table['header']]
            at_starts = {tuple(laser.table['rows'][on][header.index(c)] for c in HEADER[7:]) for on, _ in laser.table['interval_rows']}
            self.assertEqual(at_starts, {spot}, name)                   # every firing start records the same spot
            record = spot_width(laser, {})
            self.assertEqual((record['status'], record['width_um'], record['source']), ('not_computable', None, None), name)
            self.assertEqual(record['reason'], 'the spot size, type or angle changes during 1 of 3 firing intervals, '
                             'so the firing intervals have no single spot width', name)
            evidence = firing_speed(laser)
            self.assertEqual((evidence['status'], evidence['max_um_s']), ('credible', 375.), name)    # the speed is not affected
            row = distance_scale({}, [], .2, 100., evidence, record)['quantities'][0]
            self.assertEqual((row['um_max'], row['spot_widths_max']), (75., None), name)
        # The same value written again during a firing is no change. A cell of the Off record that ends a firing is not compared.
        again = rewrite(state_log(self.root/'again.csv', lines, 375, spot=('Circle', '15', '')), during_first_firing(('Circle', '15', '')))
        self.assertEqual(spot_width(load_laser(again, STATE), {})['status'], 'known')

        def in_the_off_record_that_ends_the_first_firing(rows):
            on = next(i for i, row in enumerate(rows) if row[4] == 'On')
            rows[on+1][8] = '30'
        after = rewrite(state_log(self.root/'after.csv', lines, 375, spot=('Circle', '15', '')), in_the_off_record_that_ends_the_first_firing)
        self.assertEqual(spot_width(load_laser(after, STATE), {})['status'], 'known')

    def test_declared_width_is_used_and_an_invalid_one_is_an_input_error(self):
        record = self.width('declared', ('', '15 x 15', ''), {'spot_width_along_scan_um': 5.})
        self.assertEqual(record, dict(status='known', width_um=5., source='analysis.spot_width_along_scan_um (declared)', reason=None))
        for invalid in (0, -5., True, float('nan'), float('inf'), '5'):
            with self.assertRaisesRegex(InputError, 'analysis.spot_width_along_scan_um must be a positive finite number'):
                self.width(f'bad_{invalid}', ('', '15 x 15', ''), {'spot_width_along_scan_um': invalid})


class RequestTests(Case):
    def test_both_units_in_one_source_are_a_conflict_and_one_value_on_the_command_line_replaces_the_configuration(self):
        with self.assertRaisesRegex(InputError, 'analysis.tolerance_s and analysis.tolerance_um are both given'):
            read_request({'tolerance_s': .2, 'tolerance_um': 75.})
        with self.assertRaisesRegex(InputError, '--tolerance-s and --tolerance-um are both given'):
            read_request({'tolerance_s': .2}, .2, 75.)
        with self.assertRaisesRegex(InputError, 'analysis.tolerance_s and analysis.tolerance_um are both given'):
            read_request({'tolerance_s': .2, 'tolerance_um': 75.}, .2)       # the configuration stays invalid
        self.assertEqual(read_request({'tolerance_s': .2}), dict(unit='s', value=.2, source='analysis.tolerance_s', overridden=None, saved=None))
        self.assertEqual(read_request({'tolerance_um': 75}), dict(unit='um', value=75., source='analysis.tolerance_um', overridden=None, saved=None))
        self.assertEqual(read_request({'tolerance_s': .5}, cli_um=75.),
                         dict(unit='um', value=75., source='--tolerance-um', overridden=dict(unit='s', value=.5, source='analysis.tolerance_s'), saved=None))
        self.assertEqual(read_request({'tolerance_um': 400.}, cli_s=.2),
                         dict(unit='s', value=.2, source='--tolerance-s', overridden=dict(unit='um', value=400., source='analysis.tolerance_um'), saved=None))
        self.assertEqual(read_request({}), dict(unit=None, value=None, source=None, overridden=None, saved=None))   # the estimator asks for it, as before

    def test_values_that_are_not_positive_finite_numbers_are_rejected(self):
        for value in (float('nan'), float('inf'), -float('inf'), True, False, 0, 0., -75., '75'):
            with self.assertRaisesRegex(InputError, 'analysis.tolerance_um must be a positive finite number'):
                read_request({'tolerance_um': value})
            for option in ('cli_s', 'cli_um'):
                name = '--tolerance-' + option[4:]
                with self.assertRaisesRegex(InputError, name + ' must be a positive finite number'):
                    read_request({'tolerance_s': .2}, **{option: value})
        with self.assertRaisesRegex(InputError, 'analysis.tolerance_um must be a positive finite number'):
            read_request({'tolerance_um': None})                     # a key that is present must hold a value
        self.assertEqual(estimator_settings({'tolerance_um': 75., 'spot_width_along_scan_um': 5., 'n_blocks': 10}), {'n_blocks': 10})

    def test_command_line_and_configuration_conflicts_stop_before_any_output(self):
        root = self.interval('case', tolerance_s=.2, tolerance_um=75.)
        self.assertEqual(self.error('run', root/'config.json', '--output', root/'out'),
                         'analysis.tolerance_s and analysis.tolerance_um are both given; give the tolerance in one unit')
        root = self.interval('seconds', tolerance_s=.2)
        self.assertEqual(self.error('run', root/'config.json', '--output', root/'out', '--tolerance-s', '0.2', '--tolerance-um', '75'),
                         '--tolerance-s and --tolerance-um are both given; give the tolerance in one unit')
        for option, value in (('--tolerance-s', '0'), ('--tolerance-s', '-0.2'), ('--tolerance-s', 'nan'), ('--tolerance-um', 'inf'), ('--tolerance-um', '0')):
            self.assertEqual(self.error('run', root/'config.json', '--output', root/'out', option, value), option + ' must be a positive finite number')
        for value in (0, -75., float('nan'), float('inf'), True, '75'):
            bad = self.interval(f'bad_{value}', tolerance_um=value)
            self.assertEqual(self.error('run', bad/'config.json', '--output', bad/'out'), 'analysis.tolerance_um must be a positive finite number')
            self.assertFalse((bad/'out').exists())
        self.assertFalse((root/'out').exists())
        missing = self.interval('missing')                            # no tolerance at all: the existing requirement of the estimator
        self.assertEqual(self.error('run', missing/'config.json', '--output', missing/'out', '--no-preview'),
                         'An explicit positive tolerance_s is required for the intended use')

    def test_distance_request_without_a_scan_speed_is_an_input_error_and_seconds_still_run(self):
        root = self.interval('case', tolerance_um=75.)
        message = self.error('run', root/'config.json', '--output', root/'out', '--no-preview')
        self.assertEqual(message, 'A tolerance in micrometres cannot be converted to seconds for this record: the laser log gives no scan speed '
                         'and analysis.scan_speed_um_s is not set. Give the tolerance in seconds (analysis.tolerance_s or --tolerance-s), '
                         'or declare analysis.scan_speed_um_s.')
        self.assertFalse((root/'out/report.json').exists())
        code, printed, _ = self.command('run', root/'config.json', '--output', root/'seconds', '--no-preview', '--tolerance-s', '0.2')
        self.assertEqual((code, json.loads(printed)['status']), (0, 'supported'))
        r = read(root/'seconds/report.json')
        self.assertEqual((r['distance_scale']['status'], r['distance_scale']['reason']),
                         ('not_computable', 'the laser log gives no scan speed and analysis.scan_speed_um_s is not set'))
        self.assertEqual(r['distance_scale']['quantities'][0], dict(quantity='tolerance', seconds=.2, limit_s=None, basis='time difference',
                                                                   um_min=None, um_max=None, spot_widths_min=None, spot_widths_max=None))
        self.assertEqual(r['tolerance_resolution']['overridden_configuration_value'], dict(unit='um', value=75., source='analysis.tolerance_um'))
        self.assertIn('Equivalent distances: not computable (the laser log gives no scan speed and analysis.scan_speed_um_s is not set).',
                      (root/'seconds/report.md').read_text(encoding='utf-8').splitlines())
        # With a declared speed the same request is accepted, and the declaration is its recorded source.
        declared = self.interval('declared', tolerance_um=75., scan_speed_um_s=375.)
        code, r = run(declared/'config.json', declared/'out')
        self.assertEqual((code, r['settings']['tolerance_s'], r['tolerance_resolution']['speed_source']), (0, .2, 'analysis.scan_speed_um_s (declared)'))


class RunTests(Case):
    COMPARED = ('status', 'rate_ppm', 'baseline', 'closure', 'checks', 'failed_checks', 'settings', 'alignment', 'coverage_fraction',
                'sensitivity', 'sensitivity_envelope_s', 'contradictory_high_quality_blocks', 'record_edge_lines', 'laser_log_screen',
                'record_boundary_blocks', 'duration_s', 'predicted_end_to_end_lag_change_s', 'provisional_export_eligibility')

    def test_distance_request_gives_the_result_of_its_seconds_and_is_converted_once(self):
        seconds, distance = self.scan('seconds', tolerance_s=.2), self.scan('distance', tolerance_um=75.)
        _, a = run(seconds/'config.json', seconds/'out', True)
        seen = []

        def diagnosed(ms, laser, analysis):
            seen.append(dict(analysis))
            return diagnose(ms, laser, analysis)
        with patch('lagcheck.cli.resolve_tolerance', wraps=resolve) as resolved, patch('lagcheck.cli.diagnose', side_effect=diagnosed), \
                patch('lagcheck.tolerance.firing_speed', wraps=firing_speed) as speeds:
            code, b = run(distance/'config.json', distance/'out', True)
        self.assertEqual((code, b['status'], b['correction_exported']), (0, 'supported', True))
        self.assertEqual((resolved.call_count, speeds.call_count), (1, 1))            # one conversion, on the original log
        self.assertEqual(seen, [{'alignment_offset_s': 3., 'skip_record_edge_lines': False, 'tolerance_s': .2}])   # the estimator gets seconds only
        a['provisional_export_eligibility']['tolerance_request'] = b['provisional_export_eligibility']['tolerance_request']
        for record in (a, b):
            for key in ('file',):
                record['provisional_export_eligibility']['applied_mapping'][key] = None
        for key in self.COMPARED:
            self.assertEqual(a[key], b[key], key)
        for name in ('blocks.csv', 'candidate_laser_log.csv', 'corrected_laser_log.csv', 'post_export_blocks.csv'):
            self.assertEqual(digest(seconds/'out'/name), digest(distance/'out'/name), name)
        # The read-back of the written log uses the resolved seconds; it does not derive a tolerance from the corrected log.
        self.assertEqual((b['post_export_verification']['tolerance_s'], read(distance/'out/post_export_report.json')['settings']['tolerance_s']), (.2, .2))
        record = b['tolerance_resolution']
        self.assertEqual((record['requested_unit'], record['requested_value'], record['resolved_tolerance_s'], record['speed_um_s'], record['rule']),
                         ('um', 75., .2, 375., DISTANCE_RULE))
        self.assertEqual(a['tolerance_resolution'], dict(requested_unit='s', requested_value=.2, request_source='analysis.tolerance_s',
                         overridden_configuration_value=None, resolved_tolerance_s=.2, rule=SECONDS_RULE, speed_um_s=None, speed_source=None,
                         reused_from_saved_configuration=False))
        # The saved configuration holds the seconds that were used, with the record of the request.
        saved, plain = read(distance/'out/config.json'), read(seconds/'out/config.json')
        self.assertEqual(saved['analysis'], {'alignment_offset_s': 3., 'skip_record_edge_lines': False, 'tolerance_s': .2})
        self.assertEqual(saved['tolerance_resolution'], record)
        self.assertNotIn('tolerance_resolution', plain)
        self.assertEqual(read(distance/'out/correction_manifest.json')['tolerance_resolution'], record)
        self.assertNotIn('tolerance_resolution', read(seconds/'out/correction_manifest.json'))
        text = (distance/'out/report.md').read_text(encoding='utf-8').splitlines()
        self.assertIn("Tolerance requested: 75 µm (analysis.tolerance_um). It was converted once, before the diagnosis, to 0.2 s at 375 µm s⁻¹, "
                      'the largest scan speed of the valid firing intervals. The checks and the read-back use this value in seconds.', text)
        # A rerun of the saved configuration reuses the resolved seconds and does not convert again.
        with patch('lagcheck.cli.diagnose', side_effect=diagnosed):
            code, c = run(distance/'out/config.json', distance/'again')
        self.assertEqual((code, c['rate_ppm'], c['settings']), (0, b['rate_ppm'], b['settings']))
        self.assertEqual(c['tolerance_resolution'], dict(record, reused_from_saved_configuration=True))
        self.assertEqual(seen[-1]['tolerance_s'], .2)
        changed = dict(saved, analysis=dict(saved['analysis'], tolerance_s=.3))
        (distance/'changed.json').write_text(json.dumps(changed), encoding='utf-8')
        with self.assertRaisesRegex(InputError, 'remove tolerance_resolution to give a new tolerance'):
            run(distance/'changed.json', distance/'changed')
        code, d = run(distance/'changed.json', distance/'override', tolerance_s=.3)       # an explicit new request replaces the record
        self.assertEqual((d['settings']['tolerance_s'], d['tolerance_resolution']['request_source'], d['tolerance_resolution']['reused_from_saved_configuration']),
                         (.3, '--tolerance-s', False))

    def test_command_line_value_replaces_the_tolerance_of_the_configuration_and_keeps_the_replaced_value(self):
        root = self.scan('case', tolerance_s=.5)
        code, printed, _ = self.command('run', root/'config.json', '--output', root/'um', '--no-preview', '--tolerance-um', '75')
        r = read(root/'um/report.json')
        self.assertEqual((code, r['settings']['tolerance_s']), (0, .2))
        self.assertEqual(r['tolerance_resolution'], dict(requested_unit='um', requested_value=75., request_source='--tolerance-um',
                         overridden_configuration_value=dict(unit='s', value=.5, source='analysis.tolerance_s'), resolved_tolerance_s=.2,
                         rule=DISTANCE_RULE, speed_um_s=375., speed_source=SPEED_SOURCE, reused_from_saved_configuration=False))
        self.assertEqual(read(root/'um/config.json')['analysis']['tolerance_s'], .2)
        self.assertIn('The command line replaced the tolerance of the configuration (0.5 s, analysis.tolerance_s).',
                      (root/'um/report.md').read_text(encoding='utf-8').splitlines())
        other = self.scan('other', tolerance_um=400.)
        code, r = run(other/'config.json', other/'s', tolerance_s=.2)
        self.assertEqual((code, r['settings']['tolerance_s'], r['tolerance_resolution']['requested_unit']), (0, .2, 's'))
        self.assertEqual(r['tolerance_resolution']['overridden_configuration_value'], dict(unit='um', value=400., source='analysis.tolerance_um'))
        saved = read(other/'s/config.json')
        self.assertEqual((saved['analysis'].get('tolerance_s'), 'tolerance_um' in saved['analysis']), (.2, False))

    def test_report_gives_distances_for_time_differences_only_and_spot_widths_only_from_a_declared_or_plain_width(self):
        root = self.scan('case', tolerance_s=.2, spot_width_along_scan_um=5., speeds=[225 if i % 2 else 375 for i in range(240)])
        code, r = run(root/'config.json', root/'out')
        scale, duration = r['distance_scale'], r['duration_s']
        self.assertEqual((code, scale['status'], scale['diagnostic'], scale['speed']['min_um_s'], scale['speed']['max_um_s']), (0, 'computed', 'main', 225., 375.))
        self.assertEqual(scale['spot_width'], dict(status='known', width_um=5., source='analysis.spot_width_along_scan_um (declared)', reason=None))
        rows = {row['quantity']: row for row in scale['quantities']}
        self.assertEqual(list(rows), ['tolerance', 'largest_heldout_error_single_pass', 'largest_heldout_error_last_iteration', 'sensitivity_envelope',
                                      'closure_residual_over_firing_span', 'predicted_end_to_end_lag_change', 'index_grid_lag_error_bound',
                                      'screening_half_width_over_firing_span'])
        for row in rows.values():
            self.assertAlmostEqual(row['um_min'], abs(row['seconds'])*225., places=9)
            self.assertAlmostEqual(row['um_max'], abs(row['seconds'])*375., places=9)
            self.assertAlmostEqual(row['spot_widths_max'], row['um_max']/5., places=9)
        self.assertEqual((rows['tolerance']['um_min'], rows['tolerance']['um_max'], rows['tolerance']['spot_widths_max']), (45., 75., 15.))
        self.assertEqual(rows['largest_heldout_error_single_pass']['seconds'], r['baseline']['heldout_max_s'])
        self.assertEqual(rows['sensitivity_envelope']['seconds'], r['sensitivity_envelope_s'])
        # A rate is converted only after it has been multiplied by the firing span.
        residual = rows['closure_residual_over_firing_span']
        self.assertEqual((residual['seconds'], residual['basis']), (abs(r['closure']['residual_ppm'])*1e-6*duration, 'rate difference times the firing span'))
        self.assertEqual(rows['index_grid_lag_error_bound']['limit_s'], .2/4)
        self.assertEqual(r['settings']['closure_tolerance_ppm'], .2)                 # a limit in µs/s: it has no row and is not a tolerance in seconds
        self.assertEqual(scale['not_converted'], ['coverage fractions', 'correlation peaks and margins', 'numbers of settings that gave a rate',
                                                  'check results that are true or false', 'line-change intervals flagged by the laser-log screen'])
        self.assertIn('It is a scale along a line scanned at that speed', scale['note'])
        text = (root/'out/report.md').read_text(encoding='utf-8')
        self.assertIn('- tolerance: 0.2 s; 45 to 75 µm (9 to 15 spot widths) (time difference)', text.splitlines())
        self.assertIn(f"Scan speed of the firing intervals: 225 to 375 µm s⁻¹ ({SPEED_SOURCE}; 240 firing intervals).", text.splitlines())
        self.assertNotIn('pixel', (text + json.dumps(scale)).lower())
        # The settings used for reporting never reach the estimator, and the report fields of the diagnosis are those of seconds alone.
        self.assertNotIn('spot_width_along_scan_um', r['settings'])
        self.assertEqual(read(root/'out/config.json')['analysis']['spot_width_along_scan_um'], 5.)
        plain = self.scan('plain', tolerance_s=.2, speeds=[225 if i % 2 else 375 for i in range(240)])
        _, p = run(plain/'config.json', plain/'out')
        for key in RunTests.COMPARED:
            self.assertEqual(p[key], r[key], key)
        self.assertEqual(p['distance_scale']['spot_width']['status'], 'not_computable')         # '15 x 15' without a recorded shape
        self.assertIsNone(p['distance_scale']['quantities'][0]['spot_widths_max'])
        self.assertEqual(p['distance_scale']['quantities'][0]['um_max'], 75.)

    def test_row_pitch_of_the_preview_is_not_a_spot_width(self):
        from lagcheck.preview_geometry import _pitch
        self.assertEqual((DEFAULT_PITCH_UM, _pitch([0., 0., 0.])[0]), (15., 15.))       # the display default for one line position
        laser = small_log(self.root/'log.csv', 375, header=HEADER[:7])                 # positions and speed, but no spot field
        scale = distance_scale({}, [], .2, 100., firing_speed(laser), spot_width(laser, {}))
        self.assertEqual((scale['status'], scale['spot_width']['status'], scale['quantities'][0]['um_max']), ('computed', 'not_computable', 75.))
        self.assertTrue(all(row['spot_widths_max'] is None for row in scale['quantities']))
        root = self.interval('case', tolerance_s=.2, scan_speed_um_s=375.)
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((r['distance_scale']['status'], r['distance_scale']['spot_width']['status']), ('computed', 'not_computable'))
        self.assertEqual(r['predicted_end_to_end_displacement_um'], r['predicted_end_to_end_lag_change_s']*375.)   # the existing field is as before

    def test_spot_record_runs_in_seconds_and_refuses_a_distance(self):
        root = self.scan('spots', tolerance_s=.2, move=(0., 0.))
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['status']), (0, 'supported'))
        self.assertEqual(r['distance_scale']['status'], 'not_computable')
        self.assertIn('240 of 240 firing intervals show no stage movement (spot record)', r['distance_scale']['reason'])
        message = self.error('run', root/'config.json', '--output', root/'um', '--tolerance-um', '75')
        self.assertTrue(message.startswith('A tolerance in micrometres cannot be converted to seconds for this record: 240 of 240'), message)
        self.assertTrue((root/'um/INCOMPLETE.json').exists() and (root/'um/preview_pair.png').exists())     # the raw view of an input error
        self.assertFalse((root/'um/report.json').exists())

    def test_fragments_convert_a_distance_with_their_own_speed_and_the_index_lists_seconds_and_speed(self):
        root = self.root/'fragments'
        root.mkdir()
        # One sample name in two fragments under different sequence numbers, scanned at 375 and at 225 µm/s.
        with state_log(root/'first.csv', [(10., 20.), (30., 40.), (50., 60.)], 375).open(encoding='utf-8', newline='') as f:
            rows = list(csv.reader(f))
        with state_log(root/'second.csv', [(110., 120.), (130., 140.), (150., 160.)], 225).open(encoding='utf-8', newline='') as f:
            later = list(csv.reader(f))[1:]
        later[0][0], later[0][2] = '100.000', '3'
        rows += [['90.000', 'Other', '2', '', 'Off', '0.0', '0.0', '', '15 x 15', '']] + later
        with (root/'laser.csv').open('w', encoding='utf-8', newline='') as f:
            csv.writer(f).writerows(rows)
        with (root/'signal.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['time_ms', 'element'])
            w.writerows((100*i, 1000. if (i//100) % 2 else 2.) for i in range(2000))
        config = {'ms': {'path': 'signal.csv', 'format': 'csv', 'header_row': 0, 'encoding': 'utf-8', 'delimiter': ',',
                         'time_kind': 'numeric', 'time_unit': 'ms', 'time_column': 'time_ms', 'signal_columns': ['element']},
                  'laser': dict(STATE, path='laser.csv'), 'analysis': {'tolerance_um': 75.}}
        (root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        code, r = run(root/'config.json', root/'out')
        self.assertEqual((code, r['status']), (2, 'segmented_not_jointly_verified'))
        entries = read(root/'out/segments.json')['segments']
        # The same 75 µm is 0.2 s at 375 µm/s and a third of a second at 225 µm/s. Each fragment is diagnosed on its own.
        self.assertEqual([(e['resolved_tolerance_s'], e['tolerance_conversion_speed_um_s']) for e in entries], [(.2, 375.), (75./225., 225.)])
        for entry, name in zip(entries, ('segment_001', 'segment_002')):
            record = read(root/'out'/name/'report.json')['tolerance_resolution']
            self.assertEqual((entry['resolved_tolerance_s'], entry['tolerance_conversion_speed_um_s']),
                             (record['resolved_tolerance_s'], record['speed_um_s']), name)
        self.assertEqual(r['segments'], entries)
        # A request in seconds gives every fragment those seconds, and no speed converted it.
        code, r = run(root/'config.json', root/'seconds', tolerance_s=.2)
        self.assertEqual([(e['resolved_tolerance_s'], e['tolerance_conversion_speed_um_s']) for e in r['segments']], [(.2, None)]*2)

    def test_formatted_export_reuses_the_saved_seconds(self):
        inputs = make_demo(self.root/'inputs')
        config = read(inputs/'config.json')
        del config['analysis']['tolerance_s']
        config['analysis']['tolerance_um'] = 2.
        (inputs/'config.json').write_text(json.dumps(config), encoding='utf-8')
        affine = (True, 'affine', '2026-01-01 00:00:00.000000', 'Known synthetic first stored sample')
        # The demo log sets 10 µm/s in a record above the firing starts and holds no positions. That speed is not taken
        # from the log, so the request needs the declared speed.
        with self.assertRaisesRegex(InputError, '240 of 240 firing intervals have no usable scan speed'):
            run(inputs/'config.json', self.root/'undeclared', *affine)
        config['analysis']['scan_speed_um_s'] = 10.
        (inputs/'config.json').write_text(json.dumps(config), encoding='utf-8')
        code, r = run(inputs/'config.json', self.root/'canonical', *affine)
        self.assertEqual((code, r['settings']['tolerance_s'], r['tolerance_resolution']['speed_um_s'], r['tolerance_resolution']['speed_source']),
                         (0, .2, 10., 'analysis.scan_speed_um_s (declared)'))
        with patch('lagcheck.tolerance.firing_speed', side_effect=AssertionError('no conversion after the canonical run')):
            code, result = export_iolite(self.root/'canonical', self.root/'iolite')
        self.assertEqual((code, result['post_export_verification']['tolerance_s']), (0, .2))
        self.assertEqual(read(self.root/'iolite/post_export_report.json')['settings']['tolerance_s'], .2)

    def test_distance_request_is_stated_on_the_figure_of_a_supported_result(self):
        from matplotlib.figure import Figure
        root = self.scan('case', tolerance_um=75.)
        with patch.object(Figure, 'supxlabel', autospec=True) as drawn:
            code, r = _run(root/'config.json', root/'out')
        meta = read(root/'out/preview_pair.json')
        self.assertEqual((code, r['status']), (0, 'supported'))
        self.assertEqual(meta['provisional_export_text'], ['Tolerance 0.2 s, converted before the diagnosis from 75 µm at 375 µm s⁻¹.'])
        self.assertIn(meta['provisional_export_text'][0], ' '.join(drawn.call_args_list[0].args[1].split()))


if __name__ == '__main__':
    unittest.main()
