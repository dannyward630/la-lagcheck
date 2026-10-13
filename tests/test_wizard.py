import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from lagcheck.cli import main
from lagcheck import __version__
from lagcheck.wizard import collect
from lagcheck.models import InputError
import test_cli



# Core regression tests opt out of the default rendering; visual behavior is tested separately.
from lagcheck.cli import run as _run_default_preview, main as _main_default_preview
def run(*args, **kwargs):
    if len(args) < 7:
        kwargs.setdefault('preview', False)
    elif args[6] is None:
        args = list(args); args[6] = False
    return _run_default_preview(*args, **kwargs)
def main(args=None):
    if args and args[0] == 'run' and '--preview' not in args and '--no-preview' not in args:
        args = list(args) + ['--no-preview']
    return _main_default_preview(args)

class WizardTests(unittest.TestCase):
    def test_current_version_banner_and_distinct_default_destinations(self):
        # Packaged example expectations must identify the running release too.
        expected_path = Path(__file__).resolve().parents[1] / 'examples/chromium_grating/expected_results.json'
        expected = json.loads(expected_path.read_text(encoding='utf-8'))
        for branch, section in (('automatic', 'report_subset'),
                                ('declared', 'report_subset'),
                                ('declared', 'manifest_subset')):
            with self.subTest(example_version=f'{branch}.{section}.version'):
                self.assertEqual(expected['runs'][branch][section]['version'], __version__)
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);ms,laser=self.files(root)
            destinations=[]
            for _ in range(2):
                answers=iter([str(ms),str(laser),'1','A','1','1','.2','1',''])
                messages=[]
                _,out=collect(lambda _:next(answers),messages.append)
                self.assertIn(__version__, messages[0])
                self.assertEqual(messages[-1],'Running diagnosis, protected export and complete CSV re-check; timing tolerance 0.2 s.')
                destinations.append(out)
            self.assertNotEqual(*destinations)

    def files(self,root):
        ms=root/'mass.csv';ms.write_text('example\nIntensity Vs Time,CPS\nTime [Sec],C13,Fe56\n0,1,2\n1,2,3\n2,3,4\n')
        laser=root/'laser.csv';laser.write_text('Timestamp,Laser State,Comment\n2026-01-01 00:00:00,On,A\n2026-01-01 00:00:01,Off,\n2026-01-01 00:00:02,On,B\n2026-01-01 00:00:03,Off,\n')
        return ms,laser

    def test_guided_selection_and_spatial_tolerance_conversion(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);ms,laser=self.files(root)
            answers=iter([str(ms),str(laser),'1,Fe56','B','1','2','75','375','1',str(root/'out')])
            messages=[]
            cfg,out=collect(lambda _:next(answers),messages.append)
            # The seconds shown before the run follow from the declared speed alone; the run confirms them when it reads the log.
            self.assertEqual(messages[-1],'Running diagnosis, protected export and complete CSV re-check; timing tolerance 75 µm '
                             '(0.2 s at the declared speed of 375 µm s⁻¹; the run confirms this value only after it has read the laser log).')
            self.assertEqual(cfg['ms']['signal_columns'],['C13','Fe56'])
            self.assertEqual(cfg['laser']['sample_selector'],'B')
            # The guided entry passes the request on in micrometres; the run converts it once, with the declared speed.
            self.assertEqual(cfg['analysis'],{'scan_speed_um_s':375.,'tolerance_um':75.})
            from lagcheck.laser import load_laser
            from lagcheck.tolerance import read_request,resolve
            selected=load_laser(laser,{k:v for k,v in cfg['laser'].items() if k!='path'})
            settings,record,_=resolve(read_request(cfg['analysis']),cfg['analysis'],selected)
            self.assertEqual(settings,{'tolerance_s':.2,'scan_speed_um_s':375.})
            self.assertEqual((record['requested_unit'],record['requested_value'],record['resolved_tolerance_s']),('um',75.,.2))
            self.assertEqual(cfg['export']['mode'],'rate_only')
            self.assertFalse(Path(out).exists())

    def test_distance_request_for_a_spot_record_is_refused_by_the_run_and_leaves_the_raw_view(self):
        from test_tolerance_distance import state_log
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);signal=root/'mass.csv'
            with signal.open('w',encoding='utf-8',newline='') as f:
                w=csv.writer(f)
                w.writerows([['Intensity Vs Time','CPS'],['Time [Sec]','C13']])
                w.writerows((i*.1,10 if 10<i*.1<20 else 1) for i in range(800))
            laser=state_log(root/'spots.csv',[(10.,20.),(30.,40.),(50.,60.)],375,move=(0.,0.))       # three spots: no stage movement
            out=root/'out'
            answers=iter([str(signal),str(laser),'1','1','2','s','2','75','375','1',str(out)])
            messages=[]
            cfg,destination=collect(lambda _:next(answers),messages.append)
            self.assertEqual(cfg['analysis'],{'scan_speed_um_s':375.,'tolerance_um':75.})
            self.assertIn('the run confirms this value only after it has read the laser log',messages[-1])
            printed,errors=io.StringIO(),io.StringIO()
            with patch('lagcheck.wizard.collect',return_value=(cfg,destination)),contextlib.redirect_stdout(printed),contextlib.redirect_stderr(errors):
                code=main(['wizard'])
            # A declared speed does not turn a spot record into a scan: the run refuses the distance as an input error.
            self.assertEqual((code,printed.getvalue()),(1,''))
            error=json.loads(errors.getvalue().splitlines()[-1])
            self.assertEqual(error,{'status':'input_or_tool_error','error':'A tolerance in micrometres cannot be converted to seconds for this record: '
                             '3 of 3 firing intervals show no stage movement (spot record); a time difference has no equivalent distance along a scan there. '
                             'Give the tolerance in seconds (analysis.tolerance_s or --tolerance-s).'})
            # It leaves the marker of an incomplete run with the raw view of the inputs, and neither a report nor a log.
            self.assertEqual(sorted(p.name for p in out.iterdir()),['INCOMPLETE.json','config.json','preview_pair.json','preview_pair.png'])
            self.assertEqual(json.loads((out/'INCOMPLETE.json').read_text(encoding='utf-8'))['error'],error['error'])
            view=json.loads((out/'preview_pair.json').read_text(encoding='utf-8'))
            self.assertEqual((view['view'],view['rate_fitted'],view['applied_rate_ppm']),('signal_time',False,None))

    def test_unlisted_channel_is_rejected_instead_of_summed(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);ms,laser=self.files(root)
            answers=iter([str(ms),str(laser),'Time [Sec]'])
            with self.assertRaises(InputError):collect(lambda _:next(answers),lambda _:None)

    def test_nonfinite_tolerance_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);ms,laser=self.files(root)
            answers=iter([str(ms),str(laser),'1','A','1','1','nan'])
            with self.assertRaises(InputError):collect(lambda _:next(answers),lambda _:None)

    def test_wizard_command_uses_complete_export_workflow(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);cfg=test_cli.WorkflowTests().make_case(root)
            cfg['ms']['path']=str((root/'signal.csv').resolve())
            cfg['laser']['path']=str((root/'laser.csv').resolve())
            with patch('lagcheck.wizard.collect',return_value=(cfg,str(root/'out'))),contextlib.redirect_stdout(io.StringIO()) as text:
                self.assertEqual(main(['wizard']),0)
            self.assertIn('Verified log:',text.getvalue())
            r=json.loads((root/'out/report.json').read_text(encoding='utf-8'))
            self.assertTrue(r['exported_log_diagnostic_rerun'])
            self.assertTrue(r['post_export_verification']['passed'])
            self.assertTrue(Path(r['usable_file']).is_file())


if __name__=='__main__':unittest.main()
