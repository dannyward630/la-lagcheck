import contextlib
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
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);ms,laser=self.files(root)
            destinations=[]
            for _ in range(2):
                answers=iter([str(ms),str(laser),'1','A','1','1','.2','1',''])
                messages=[]
                _,out=collect(lambda _:next(answers),messages.append)
                self.assertIn(__version__, messages[0])
                destinations.append(out)
            self.assertNotEqual(*destinations)

    def files(self,root):
        ms=root/'mass.csv';ms.write_text('example\nIntensity Vs Time,CPS\nTime [Sec],C13,Fe56\n0,1,2\n1,2,3\n2,3,4\n')
        laser=root/'laser.csv';laser.write_text('Timestamp,Laser State,Comment\n2026-01-01 00:00:00,On,A\n2026-01-01 00:00:01,Off,\n2026-01-01 00:00:02,On,B\n2026-01-01 00:00:03,Off,\n')
        return ms,laser

    def test_guided_selection_and_spatial_tolerance_conversion(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);ms,laser=self.files(root)
            answers=iter([str(ms),str(laser),'1,Fe56','B','1','2','75','375','1',str(root/'out')])
            cfg,out=collect(lambda _:next(answers),lambda _:None)
            self.assertEqual(cfg['ms']['signal_columns'],['C13','Fe56'])
            self.assertEqual(cfg['laser']['sample_selector'],'B')
            self.assertEqual(cfg['analysis'],{'tolerance_s':.2,'scan_speed_um_s':375.})
            self.assertEqual(cfg['export']['mode'],'rate_only')
            self.assertFalse(Path(out).exists())

    def test_unlisted_channel_is_rejected_instead_of_summed(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);ms,laser=self.files(root)
            answers=iter([str(ms),str(laser),'Time [Sec]'])
            with self.assertRaises(InputError):collect(lambda _:next(answers),lambda _:None)

    def test_nonfinite_tolerance_rejected(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);ms,laser=self.files(root)
            answers=iter([str(ms),str(laser),'1','A','1','1','nan'])
            with self.assertRaises(InputError):collect(lambda _:next(answers),lambda _:None)

    def test_wizard_command_uses_complete_export_workflow(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
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
