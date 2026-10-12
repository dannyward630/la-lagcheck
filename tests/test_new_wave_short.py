import csv
import json
from pathlib import Path
import numpy as np
import pytest
from lagcheck.models import InputError, LaserData, SignalData
from lagcheck.laser import load_laser, export_corrected
from lagcheck.processing import inspect_segments
from lagcheck.ms import load_signal
from lagcheck.estimator import diagnose
from lagcheck.cli import run


HEADER = ['Timestamp','Sequence Number','SubPoint Number','Vertix Number','Comment',
          'X(um)','Y(um)','Scan Velocity (um/s)','Laser State']
CFG = dict(format='new_wave_laser_log', time_kind='numeric', time_unit='s', sample_scope='all')


def log(tmp_path, rows):
    path=tmp_path/'laser.csv'
    with path.open('w', newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(HEADER);w.writerows(rows)
    return path


def records(n=3, spots=False):
    rows=[]
    for i in range(n):
        start=10*i
        rows += [[start,str(i+1),'1','', 'sample' if i%2==0 else 'standard', i*10,0,'','Off'],
                 [start+1,'','','','',i*10,0,'','On'],
                 [start+1.1,'','','2','',i*10,0,'1','On'],
                 [start+5,'','','','',i*10 if spots else i*10+4,0,'','Off']]
    return rows


def test_new_wave_preserves_states_sequences_and_source_records(tmp_path):
    path=log(tmp_path,records())
    laser=load_laser(path,CFG)
    assert len(laser.lines)==3
    assert laser.metadata['new_wave']['state_update_count']==3
    assert laser.metadata['acquisition_identity']['sequence_ids']==['1','2','3']
    assert [x['start_record'] for x in laser.metadata['new_wave']['intervals']]==[2,6,10]
    assert len(inspect_segments(path,CFG)['segments'])==1
    assert laser.lines[0]['Start Time']==1
    assert laser.lines[0]['End Time']==5


def test_new_wave_disjoint_named_selection_does_not_bridge_standard(tmp_path):
    path=log(tmp_path,records())
    cfg=dict(CFG);cfg.pop('sample_scope');cfg['sample_selector']='sample'
    laser=load_laser(path,cfg)
    assert [(v['Start Time'],v['End Time']) for v in laser.lines]==[(1,5),(21,25)]
    assert laser.metadata['input_processing']['selected_record_numbers']==[1,2,3,4,9,10,11,12]
    assert len(inspect_segments(path,cfg)['segments'])==1
    receipt=export_corrected(laser,tmp_path/'corrected.csv',150)
    parsed=load_laser(tmp_path/'corrected.csv',receipt['output_config'])
    assert len(parsed.lines)==2
    assert receipt['non_time_cells_verified'] and receipt['firing_intervals_verified']


@pytest.mark.parametrize('change', ['sequence','sample','decreasing_time','unclosed'])
def test_new_wave_rejects_ambiguous_firing_changes(tmp_path,change):
    rows=records(1)
    if change=='sequence': rows[2][1]='2'
    elif change=='sample': rows[2][4]='other'
    elif change=='decreasing_time': rows[2][0]=.5
    else: rows=rows[:-1]
    with pytest.raises(InputError): load_laser(log(tmp_path,rows),CFG)


def test_new_wave_requires_explicit_all_sample_scope(tmp_path):
    cfg=dict(CFG);cfg.pop('sample_scope')
    with pytest.raises(InputError,match='Multiple samples'): load_laser(log(tmp_path,records()),cfg)


def test_generic_state_reader_still_rejects_repeated_on(tmp_path):
    cfg=dict(CFG,format='laser_state_csv');cfg.pop('sample_scope')
    with pytest.raises(InputError,match='Consecutive On'): load_laser(log(tmp_path,records(1)),cfg)


def ms_file(tmp_path,trailer):
    path=tmp_path/'ms.csv'
    path.write_text('Intensity Vs Time,Counts\nTime [Sec],Ca43\n0,1\n1,2\n2,3\n'+trailer,encoding='utf-8')
    return path


def test_masshunter_delimited_footer_padding_is_explicitly_recorded(tmp_path):
    ms=load_signal(ms_file(tmp_path,',\n,\nPrinted: 2011-06-15,\n'),
                   dict(format='agilent_masshunter_csv',signal_columns=['Ca43']))
    assert ms.metadata['recognized_footer_padding_physical_rows']==[6,7]
    assert ms.metadata['row_count']==3
    assert ms.signal.tolist()==[1,2,3]


@pytest.mark.parametrize('trailer',[',\n',',\n3,4\n',',\nPrinted: something,nonblank\n'])
def test_masshunter_does_not_discard_missing_rows_as_padding(tmp_path,trailer):
    with pytest.raises(InputError):
        load_signal(ms_file(tmp_path,trailer),dict(format='agilent_masshunter_csv',signal_columns=['Ca43']))


def short_data(n=72):
    rng=np.random.default_rng(10)
    lines=[];cursor=5.
    for _ in range(n):
        duration=rng.uniform(5,18)
        lines.append({'Start Time':cursor,'End Time':cursor+duration})
        cursor+=duration+rng.uniform(3,10)
    t=np.arange(0,cursor+30,.1);y=np.ones(len(t))*2
    for row in lines:
        a=5+1.0001*(row['Start Time']-5)+3
        b=5+1.0001*(row['End Time']-5)+3
        y+=1000*np.maximum(0,np.minimum(t+.05,b)-np.maximum(t-.05,a))/.1
    return SignalData(t,y),LaserData(lines,metadata={'format':'new_wave_laser_log'})


def test_short_record_keeps_gates_and_reports_unexecuted_sensitivity():
    ms,laser=short_data()
    r,blocks=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':3.})
    assert r['status']=='short_acquisition_conditional'
    assert abs(r['rate_ppm']-100)<20
    assert len(blocks)==6
    short=r['short_acquisition']
    assert short['minimum_intervals_per_block']==12 and short['peak_min']==.25 and short['margin_min']==.01
    assert short['unavailable_sensitivity_configs']==['n_blocks=8','n_blocks=12','n_blocks=15']
    assert r['standard_diagnostic']['reason']=='insufficient_laser_intervals'
    assert short['check_execution']['baseline_heldout']['state']=='pass'
    assert short['check_execution']['all_sensitivity_configs_estimable']['state']=='incomplete'


def test_strict_switch_and_ordinary_input_preserve_standard_behavior():
    ms,laser=short_data()
    settings={'tolerance_s':.2,'alignment_offset_s':3.,'short_acquisition':False}
    r,_=diagnose(ms,laser,settings)
    assert r['reason']=='insufficient_laser_intervals' and 'short_acquisition' not in r
    laser.metadata={}
    r,_=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':3.})
    assert r['reason']=='insufficient_laser_intervals' and 'short_acquisition' not in r


def test_fewer_than_four_blocks_does_not_manufacture_rate():
    ms,laser=short_data(36)
    r,_=diagnose(ms,laser,{'tolerance_s':.2,'alignment_offset_s':3.})
    assert r['rate_ppm'] is None
    assert r['short_acquisition']['executed_blocks']==0
    execution=r['short_acquisition']['check_execution']
    # The laser-log screen needs no fit; every other check stays unexecuted.
    assert all(v['state']=='not_executed' for k,v in execution.items() if k!='no_flagged_line_change_interval')
    assert execution['no_flagged_line_change_interval']['source']=='laser_log_screen'


def test_short_unidentified_peak_remains_failed_during_conditional_fit(monkeypatch):
    from lagcheck import estimator
    ms,laser=short_data()
    monkeypatch.setattr(estimator,'_coarse',lambda *args: dict(ok=False,offset_s=3.,peak=.9,
        margin=.009,boundary=False,reason='coarse_alignment_check'))
    r,blocks=diagnose(ms,laser,{'tolerance_s':.2})
    assert len(blocks)==6 and r['rate_ppm'] is not None
    assert r['alignment']['ok'] is False
    assert r['alignment']['source']=='unverified_automatic_peak_candidate'
    assert r['checks']['initial_alignment_identifiable'] is False
    assert r['short_acquisition']['check_execution']['initial_alignment_identifiable']['state']=='fail'
    assert r['status']=='short_acquisition_conditional'


def test_short_default_point_views_even_without_rate(tmp_path,monkeypatch):
    from lagcheck import cli
    real_diagnose=cli.diagnose
    def conditional(*args):
        r,b=real_diagnose(*args)
        r['alignment']=dict(ok=False,offset_s=3.,source='unverified_automatic_peak_candidate')
        return r,b
    monkeypatch.setattr(cli,'diagnose',conditional)
    rows=records(6,spots=True)
    # Actual New Wave files can repeat the timestamp with stale X/Y metadata.
    # Those updates must not create fictitious trajectory vertices or drop a spot.
    rows[2][0]=rows[1][0]
    rows[2][5]=99999
    laser=log(tmp_path,rows)
    ms=tmp_path/'signal.csv'
    ms.write_text('t,Fe56\n'+'\n'.join(f'{i/10},{1+(i//10)%10}' for i in range(650)),encoding='utf-8')
    config={'ms':{'path':str(ms),'format':'csv','header_row':0,'delimiter':',','encoding':'utf-8','time_column':'t','time_kind':'numeric','time_unit':'s','signal_columns':['Fe56']},
            'laser':dict(CFG,path=str(laser)),'analysis':{'tolerance_s':.2,'alignment_offset_s':0.}}
    code,report=run(config,tmp_path/'out',export=True)
    assert code==2 and report['rate_ppm'] is None and not report['correction_exported']
    out=tmp_path/'out'
    meta=json.loads((out/'new_wave_views.json').read_text())
    assert meta['view']=='spot_position_map' and meta['interval_count']==6
    assert meta['geometry_errors']==[] and meta['intermediate_on_positions_used'] is False
    assert meta['panels']==['Unverified automatic alignment; no rate correction']
    preview=json.loads((out/'preview_pair.json').read_text())
    assert preview['alignment_offset_s']==3.
    assert preview['panels']==['Unverified automatic alignment; no rate correction']
    assert preview['scan_speed_um_s'] is None and 'positioning' in preview['scan_speed_source']
    assert preview['span_displacement_um'] is None
    assert len(meta['panels'])==1 and meta['spatial_accuracy_verified'] is False
    for name in ['preview_pair.png','new_wave_spatial.png','new_wave_signal_time.png','new_wave_block_lags.png']:
        assert (out/name).stat().st_size>1000
