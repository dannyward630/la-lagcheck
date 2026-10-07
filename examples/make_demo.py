"""Create a reproducible 100 µs s⁻¹ example without touching existing files."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def make_demo(folder):
    folder=Path(folder)
    folder.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(148)
    intervals=[]
    cursor=5.
    for i in range(240):
        length=float(rng.uniform(5,18))
        intervals.append((cursor,cursor+length))
        cursor+=length+float(rng.uniform(3,10))
    time=np.arange(0,cursor+30,.1)
    signal=np.full(len(time),2.)
    for start,end in intervals:
        start=5+(start-5)*1.0001+3
        end=5+(end-5)*1.0001+3
        signal+=1000*np.maximum(0,np.minimum(time+.05,end)-np.maximum(time-.05,start))/.1
    signal+=rng.normal(0,.05,len(time))
    with (folder/'signal.csv').open('x',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['time_s','signal']);w.writerows(zip(time,signal))
    with (folder/'laser.csv').open('x',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['start_s','end_s']);w.writerows(intervals)
    config={'ms':{'path':'signal.csv','format':'csv','header_row':0,'encoding':'utf-8',
                  'delimiter':',','time_kind':'numeric','time_unit':'s','time_column':'time_s',
                  'signal_columns':['signal']},
            'laser':{'path':'laser.csv','format':'interval_csv','time_kind':'numeric','time_unit':'s',
                     'start_column':'start_s','end_column':'end_s'},
            'analysis':{'tolerance_s':.2,'alignment_offset_s':8.}}
    # Independent rebasing: MS starts at 0 s; laser starts at 5 s. Physical lag is 3 s.
    (folder/'config.json').write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    (folder/'truth.json').write_text(json.dumps({'rate_ppm':100.,'physical_offset_s':3.,
        'rebased_alignment_offset_s':8.,'n_intervals':240,'seed':148},indent=2)+'\n',encoding='utf-8')
    print(folder.resolve())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True)
    make_demo(p.parse_args().output)
