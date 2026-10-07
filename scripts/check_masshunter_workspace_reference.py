"""Read-only check of a CSV time column against an existing iolite workspace.

Requires h5py for this local evidence check; it is not part of the correction
tool's runtime. An observed epoch/calendar difference is reported, never applied
as an inferred universal timezone or clock correction.
"""
import argparse
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re
import h5py
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--csv',required=True);p.add_argument('--workspace',required=True)
    p.add_argument('--time-dataset',default='Inputs/Fe56_time')
    p.add_argument('--output',required=True)
    a=p.parse_args();signal=Path(a.csv);workspace=Path(a.workspace)
    before={'csv':sha(signal),'workspace':sha(workspace)}
    with signal.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.reader(f))
    at=next(i for i,row in enumerate(rows) if row and row[0].strip()=='Time [Sec]')
    acquired=next(row[0] for row in rows[:at] if row and row[0].strip().startswith('Acquired'))
    match=re.search(r'Acquired\s*:\s*(.*?)\s+using',acquired)
    header=datetime.strptime(match.group(1),'%d-%b-%y %H:%M:%S')
    times=np.array([float(row[0]) for row in rows[at+1:] if len(row)==len(rows[at]) and row[0].strip()])
    with h5py.File(workspace,'r') as f:stored=f[a.time_dataset][:]
    if len(stored)<len(times):raise ValueError('Workspace has fewer points than CSV; no direct prefix correspondence')
    offsets=stored[:len(times)]-times
    local_header_seconds=(header-datetime(1970,1,1)).total_seconds()
    first_local=header+timedelta(seconds=float(times[0]))
    result=dict(csv=str(signal.resolve()),workspace=str(workspace.resolve()),time_dataset=a.time_dataset,
                input_sha256=before,csv_points=len(times),workspace_points=len(stored),
                extra_workspace_points=len(stored)-len(times),header_time_naive=header.isoformat(' '),
                first_csv_time_s=float(times[0]),first_sample_local_naive=first_local.isoformat(' ',timespec='microseconds'),
                workspace_first_epoch_s=float(stored[0]),
                workspace_first_utc=datetime.fromtimestamp(float(stored[0]),timezone.utc).isoformat(timespec='microseconds'),
                workspace_minus_csv_time_min_s=float(offsets.min()),workspace_minus_csv_time_max_s=float(offsets.max()),
                offset_range_s=float(np.ptp(offsets)),
                recorded_epoch_minus_naive_header_s=float(offsets[0]-local_header_seconds),
                interpretation='Evidence for this file pair only. The calendar/epoch difference is not an extra timing correction and does not establish a calibrated clock.',
                input_unchanged=before=={'csv':sha(signal),'workspace':sha(workspace)})
    if not result['input_unchanged']:raise ValueError('Input changed during reference check')
    with Path(a.output).open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
