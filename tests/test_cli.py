import csv
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from lagcheck.cli import run,main,digest
from lagcheck.ms import load_signal
from lagcheck.laser import load_laser
from lagcheck.estimator import diagnose
from test_estimator import synthetic



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

class WorkflowTests(unittest.TestCase):
    def make_case(self,root,step=0.):
        ms,laser=synthetic(100.,3.,step)
        with (root/'signal.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['time_ms','element','not_signal'])
            w.writerows(zip(ms.time*1000,ms.signal,np.ones(len(ms.time))*999))
        with (root/'laser.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['start','end','x'])
            w.writerows((l['Start Time'],l['End Time'],i) for i,l in enumerate(laser.lines))
        config={'ms':{'path':'signal.csv','format':'csv','header_row':0,'encoding':'utf-8',
                      'delimiter':',','time_kind':'numeric','time_unit':'ms','time_column':'time_ms',
                      'signal_columns':['element']},
                'laser':{'path':'laser.csv','format':'interval_csv','header_row':0,'encoding':'utf-8',
                         'delimiter':',','time_kind':'numeric','time_unit':'s',
                         'start_column':'start','end_column':'end'},
                'analysis':{'tolerance_s':.2,'alignment_offset_s':8.}}
        # Importers rebase each file at its first timestamp: laser origin is 5 s, MS origin 0 s.
        # Hence the 3 s physical offset is represented by 8 s between the rebased coordinates.
        (root/'config.json').write_text(json.dumps(config),encoding='utf-8')
        return config

    def test_end_to_end_export_and_remeasurement(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);config=self.make_case(root)
            before={n:digest(root/n) for n in ('signal.csv','laser.csv')}
            code,r=run(root/'config.json',root/'out',True)
            self.assertEqual(code,0,r)
            self.assertTrue(r['correction_exported'])
            self.assertFalse(r['importer_roundtrip_tested'])
            self.assertTrue(r['exported_log_diagnostic_rerun'])
            self.assertTrue(r['post_export_verification']['passed'])
            self.assertFalse(r['post_export_verification']['automatic_realignment_used'])
            self.assertTrue((root/'out/corrected_laser_log.csv').is_file())
            self.assertAlmostEqual(r['rate_ppm'],100.,delta=3.)
            export_config=r['export_verification'].get('output_config',config['laser'])
            export_config={k:v for k,v in export_config.items() if k!='path'}
            corrected=load_laser(root/'out/corrected_laser_log.csv',export_config)
            ms=load_signal(root/'signal.csv',{k:v for k,v in config['ms'].items() if k!='path'})
            rr,_=diagnose(ms,corrected,config['analysis'])
            self.assertLess(abs(rr['rate_ppm']),3.)
            for n,h in before.items(): self.assertEqual(digest(root/n),h)
            self.assertEqual(main(['run',str(root/'config.json'),'--output',str(root/'out'),'--export']),1)

    def test_qc_failure_does_not_export(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);self.make_case(root,step=2.)
            code,r=run(root/'config.json',root/'out',True)
            self.assertEqual(code,2)
            self.assertFalse(r['correction_exported'])
            self.assertFalse((root/'out/corrected_laser_log.csv').exists())
            self.assertTrue((root/'out/report.json').exists())

    def test_explicit_netcdf_mapping_and_cli_inspection(self):
        from scipy.io import netcdf_file
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder);cfg=self.make_case(root)
            values=np.loadtxt(root/'signal.csv',delimiter=',',skiprows=1)
            with netcdf_file(root/'signal.nc','w') as nc:
                nc.createDimension('sample',len(values))
                for col,name in enumerate(('clock','element','unused')):
                    v=nc.createVariable(name,'d',('sample',))
                    v[:]=values[:,col]
                    v.units='ms' if name=='clock' else 'counts/s'
            cfg['ms']={'path':'signal.nc','format':'netcdf','time_variable':'clock',
                       'time_unit':'ms','signal_variables':['element']}
            (root/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
            before=digest(root/'signal.nc')
            buffer=io.StringIO()
            with contextlib.redirect_stdout(buffer):
                self.assertEqual(main(['inspect-netcdf',str(root/'signal.nc')]),0)
            self.assertIn('clock',json.loads(buffer.getvalue())['variables'])
            code,r=run(root/'config.json',root/'out',True)
            self.assertEqual(code,0,r)
            self.assertAlmostEqual(r['rate_ppm'],100.,delta=3.)
            self.assertTrue(r['correction_exported'])
            self.assertEqual(digest(root/'signal.nc'),before)

    def test_invalid_file_is_tool_input_error(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as folder:
            root=Path(folder)
            (root/'config.json').write_text('{"ms": {}, "laser": {}, "analysis": {}}')
            self.assertEqual(main(['run',str(root/'config.json'),'--output',str(root/'out')]),1)
            self.assertFalse((root/'out').exists())

if __name__=='__main__':unittest.main()
