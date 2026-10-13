"""Preview pair: drawn after report.json, never changing the diagnosis, the report or the export."""
import contextlib
import csv
from datetime import datetime, timedelta
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from lagcheck.cli import run, main, digest
from lagcheck.models import InputError
from lagcheck.preview import render_preview_pair, rasterize_true_y
import test_cli
from test_estimator import synthetic

try:
    import matplotlib  # noqa: F401
    HAVE_MATPLOTLIB = True
except ImportError:
    HAVE_MATPLOTLIB = False


def file_hashes(folder):
    return {p.name: digest(p) for p in Path(folder).iterdir() if p.is_file()}


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


@unittest.skipUnless(HAVE_MATPLOTLIB, 'matplotlib (optional extra "preview") is not installed')
class PreviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def case(self, **analysis):
        """The interval-CSV fixture of test_cli: 240 lines, 100 µs/s, 8 s rebased alignment."""
        config = test_cli.WorkflowTests().make_case(self.root)
        config['analysis'].update(analysis)
        (self.root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        return config

    def esl_case(self):
        """An ESL-style state log with X(um)/Y(um) positions and an inherited scan velocity."""
        ms, laser = synthetic(100., 3.)
        with (self.root/'signal.csv').open('w', newline='') as f:
            w = csv.writer(f); w.writerow(['time_ms', 'element'])
            w.writerows(zip(ms.time*1000, ms.signal))
        base = datetime(2026, 1, 1)
        stamp = lambda s: (base+timedelta(seconds=s)).strftime('%Y-%m-%d %H:%M:%S.%f')  # noqa: E731
        with (self.root/'laser.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['Timestamp', 'Comment', 'Sequence Number', 'Scan Velocity (um/s)', 'Laser State', 'X(um)', 'Y(um)'])
            w.writerow([stamp(0), 'Synthetic A', '1', '50', 'Off', '0', '0'])
            for i, line in enumerate(laser.lines):
                length = 50.*(line['End Time']-line['Start Time'])
                xs, xe = (0., length) if i % 2 == 0 else (length, 0.)   # serpentine raster
                w.writerow([stamp(line['Start Time']), '', '', '', 'On', f'{xs:.1f}', f'{20.*i:.1f}'])
                w.writerow([stamp(line['End Time']), '', '', '', 'Off', f'{xe:.1f}', f'{20.*i:.1f}'])
        config = {'ms': {'path': 'signal.csv', 'format': 'csv', 'header_row': 0, 'encoding': 'utf-8',
                         'delimiter': ',', 'time_kind': 'numeric', 'time_unit': 'ms', 'time_column': 'time_ms',
                         'signal_columns': ['element']},
                  'laser': {'path': 'laser.csv', 'format': 'laser_state_csv', 'time_kind': 'datetime',
                            'sample_selector': 'Synthetic A'},
                  'analysis': {'tolerance_s': .2, 'alignment_offset_s': 3.}}
        (self.root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        return config

    def test_fitted_rate_gives_both_panels_and_finite_numbers(self):
        self.case(scan_speed_um_s=50.)
        code, r = run(self.root/'config.json', self.root/'out', preview=True)
        self.assertIn(code, (0, 2))
        self.assertEqual(r['preview']['status'], 'written', r['preview'])
        self.assertGreater((self.root/'out/preview_pair.png').stat().st_size, 1000)
        meta = read_json(self.root/'out/preview_pair.json')
        self.assertEqual(meta['panels'], ['one constant', 'correction: acquisition_estimate'])
        self.assertTrue(meta['rate_fitted'])
        self.assertEqual(meta['channel'], 'element')
        self.assertEqual(meta['channel_note'], 'same channel as the diagnosis')
        for key in ('constant_lag_s', 'baseline_rate_ppm', 'first_line_lag_s', 'last_line_lag_s', 'lag_span_s',
                    'scan_speed_um_s', 'first_line_displacement_um', 'last_line_displacement_um',
                    'span_displacement_um', 'predicted_first_to_last_displacement_um', 'pixel_pitch'):
            self.assertIsInstance(meta[key], float, key)
            self.assertTrue(np.isfinite(meta[key]), key)
        self.assertAlmostEqual(meta['baseline_rate_ppm'], 100., delta=3.)
        # Rebased coordinates: laser origin 5 s, MS origin 0 s, physical lag 3 s, so the first line lags 8 s.
        self.assertAlmostEqual(meta['first_line_lag_s'], 8., delta=.1)
        self.assertAlmostEqual(meta['lag_span_s'], meta['applied_rate_ppm']*1e-6*meta['first_to_last_line_start_span_s'], places=6)
        self.assertAlmostEqual(meta['predicted_first_to_last_displacement_um'],
                               50.*meta['applied_rate_ppm']*1e-6*meta['first_to_last_line_start_span_s'], places=6)
        self.assertAlmostEqual(meta['span_displacement_um'], meta['predicted_first_to_last_displacement_um'], places=6)
        low, high = sorted((meta['first_line_lag_s'], meta['last_line_lag_s']))
        self.assertTrue(low <= meta['constant_lag_s'] <= high)
        self.assertEqual(meta['scan_speed_source'], 'analysis.scan_speed_um_s')
        self.assertEqual(meta['grid_shape'][0], 240)
        self.assertAlmostEqual(meta['pixel_pitch'], .1)
        self.assertEqual(meta['pixel_pitch_unit'], 's')
        self.assertEqual(meta['colour_scale']['log10_percentiles'], [2, 98])
        self.assertLess(meta['colour_scale']['vmin'], meta['colour_scale']['vmax'])
        self.assertTrue(meta['input_unchanged'])

    def test_without_scan_speed_rows_are_time_based_and_displacements_absent(self):
        self.case()
        code, r = run(self.root/'config.json', self.root/'out', preview=True)
        meta = read_json(self.root/'out/preview_pair.json')
        self.assertTrue(meta['rate_fitted'])
        self.assertIsNone(meta['scan_speed_um_s'])
        self.assertIsNone(meta['predicted_first_to_last_displacement_um'])
        self.assertEqual(meta['pixel_pitch_unit'], 's')
        self.assertAlmostEqual(meta['pixel_pitch'], .1, places=6)
        self.assertEqual(meta['grid_shape'][0], 240)

    def test_summed_channels_default_to_first_stored_channel_and_say_so(self):
        config = self.case()
        config['ms']['signal_columns'] = ['element', 'not_signal']
        (self.root/'config.json').write_text(json.dumps(config), encoding='utf-8')
        code, r = run(self.root/'config.json', self.root/'out', preview=True)
        meta = read_json(self.root/'out/preview_pair.json')
        self.assertEqual(meta['channel'], 'element')
        self.assertEqual(meta['stored_channels'], ['element', 'not_signal'])
        self.assertIn('first stored channel', meta['channel_note'])
        # An explicitly chosen flat channel (constant 999) must still draw, with a widened scale.
        summary = render_preview_pair(self.root/'out', 'not_signal', name='flat')
        self.assertEqual(summary['channel'], 'not_signal')
        self.assertIn('the sum of', summary['channel_note'])
        self.assertTrue((self.root/'out/flat.png').is_file())
        self.assertLess(summary['colour_scale']['vmin'], summary['colour_scale']['vmax'])

    def test_no_fitted_rate_gives_left_panel_only(self):
        self.case()
        with (self.root/'laser.csv').open(newline='') as f:
            rows = list(csv.reader(f))
        with (self.root/'laser.csv').open('w', newline='') as f:
            csv.writer(f).writerows(rows[:31])      # 30 intervals: diagnosis stops before any fit
        code, r = run(self.root/'config.json', self.root/'out', preview=True)
        self.assertEqual(code, 2)
        report = read_json(self.root/'out/report.json')
        self.assertEqual(report['reason'], 'insufficient_laser_intervals')
        self.assertNotIn('baseline', report)
        self.assertEqual(r['preview']['status'], 'written', r['preview'])
        meta = read_json(self.root/'out/preview_pair.json')
        self.assertFalse(meta['rate_fitted'])
        self.assertEqual(meta['panels'], ['original / declared alignment'])
        self.assertEqual(meta['constant_lag_s'], 8.)
        self.assertEqual(meta['constant_source'], 'explicit_user_offset')
        self.assertIsNone(meta['baseline_rate_ppm'])
        self.assertIsNone(meta['first_line_lag_s'])
        self.assertEqual(meta['grid_shape'][0], 30)
        self.assertGreater((self.root/'out/preview_pair.png').stat().st_size, 1000)

    def test_preview_does_not_alter_report_or_other_outputs(self):
        self.case()
        code_a, _ = run(self.root/'config.json', self.root/'plain', preview=False)
        code_b, r = run(self.root/'config.json', self.root/'with_preview', preview=True)
        self.assertEqual(code_a, code_b)
        plain, drawn = file_hashes(self.root/'plain'), file_hashes(self.root/'with_preview')
        self.assertEqual(set(drawn)-set(plain), {'preview_pair.png', 'preview_pair.json'})
        self.assertEqual({k: drawn[k] for k in plain}, plain)   # report.json, blocks, config: identical bytes
        self.assertNotIn('preview', read_json(self.root/'with_preview/report.json'))
        # Drawing later into the plain run changes none of its existing files either.
        before = file_hashes(self.root/'plain')
        render_preview_pair(self.root/'plain')
        after = file_hashes(self.root/'plain')
        self.assertEqual({k: after[k] for k in before}, before)
        self.assertFalse((self.root/'plain/INCOMPLETE.json').exists())
        with self.assertRaisesRegex(InputError, 'already exists'):
            render_preview_pair(self.root/'plain')
        self.assertEqual(file_hashes(self.root/'plain'), after)

    def test_preview_failure_leaves_the_completed_run_intact(self):
        self.case()
        errors = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()) as console, contextlib.redirect_stderr(errors):
            code = main(['run', str(self.root/'config.json'), '--output', str(self.root/'out'), '--preview', 'no_such_column'])
        self.assertEqual(code, 1)
        summary = json.loads(console.getvalue())
        self.assertEqual(summary['preview']['status'], 'fallback_written')
        self.assertIn('no_such_column', summary['preview']['error'])
        self.assertIn('preview', errors.getvalue())
        self.assertTrue((self.root/'out/report.json').is_file())
        self.assertFalse((self.root/'out/INCOMPLETE.json').exists())
        self.assertTrue((self.root/'out/preview_pair.png').exists())
        self.assertEqual(read_json(self.root/'out/preview_pair.json')['view'],'signal_time')

    def test_cli_flag_with_channel_and_preview_subcommand(self):
        self.case()
        with contextlib.redirect_stdout(io.StringIO()) as console:
            code = main(['run', str(self.root/'config.json'), '--output', str(self.root/'out'), '--preview', 'element'])
        self.assertIn(code, (0, 2))
        self.assertEqual(json.loads(console.getvalue())['preview']['status'], 'written')
        self.assertTrue((self.root/'out/preview_pair.png').is_file())
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(['preview', str(self.root/'out')]), 1)    # exists: never overwritten
        with contextlib.redirect_stdout(io.StringIO()) as console:
            self.assertEqual(main(['preview', str(self.root/'out'), '--channel', 'element', '--name', 'second']), 0)
        self.assertEqual(json.loads(console.getvalue())['channel'], 'element')
        self.assertTrue((self.root/'out/second.png').is_file())
        self.assertTrue((self.root/'out/second.json').is_file())

    def test_esl_log_positions_define_the_grid(self):
        self.esl_case()
        code, r = run(self.root/'config.json', self.root/'out', preview=True)
        self.assertIn(code, (0, 2))
        self.assertEqual(r['preview']['status'], 'written', r['preview'])
        meta = read_json(self.root/'out/preview_pair.json')
        self.assertTrue(meta['rate_fitted'])
        self.assertTrue(meta['geometry_source'].startswith('laser log columns'))
        self.assertFalse(meta['axes_swapped'])
        self.assertEqual(meta['pixel_pitch'], 20.)
        self.assertEqual(meta['grid_shape'][0], 240)
        self.assertEqual(meta['scan_speed_um_s'], 50.)
        self.assertIn('Scan Velocity', meta['scan_speed_source'])
        self.assertTrue(meta['scan_speed_uniform_in_log'])
        self.assertAlmostEqual(meta['first_line_lag_s'], 3., delta=.1)
        self.assertAlmostEqual(meta['baseline_rate_ppm'], 100., delta=3.)
        self.assertAlmostEqual(meta['predicted_first_to_last_displacement_um'],
                               50.*meta['applied_rate_ppm']*1e-6*meta['first_to_last_line_start_span_s'], places=6)
        # The longest line is 18 s at 50 µm/s = 900 µm, so the grid has about 46 columns of 20 µm.
        self.assertGreater(meta['grid_shape'][1], 30)
        self.assertLess(meta['grid_shape'][1], 60)

    def test_rasterizer_places_a_known_line(self):
        t = np.arange(0., 10., .1)
        chan = np.where((t >= 2.) & (t < 4.), 100., 1.)
        L = [{'Start Time': 1., 'End Time': 3., 'StartX': 0., 'EndX': 200., 'StartY': 0.}]
        img = rasterize_true_y(t, chan, L, np.array([1.]), 0., 21, 10., 0., 1, 10.)
        self.assertEqual(img.shape, (1, 21))
        self.assertTrue(np.all(img[0, :20] > 50.))
        self.assertTrue(np.isnan(img[0, 20]))
        shifted = rasterize_true_y(t, chan, L, np.array([0.]), 0., 21, 10., 0., 1, 10.)
        self.assertTrue(np.all(shifted[0, :10] < 2.))
        self.assertTrue(np.all(shifted[0, 10:20] > 50.))


if __name__ == '__main__':
    unittest.main()
