"""A small guided entry for the tested MassHunter and state-log pathway."""
import csv
from datetime import datetime, timezone
import math
from pathlib import Path
import re
from uuid import uuid4

from . import __version__
from .models import InputError


def _path(value):
    path = Path(value.strip().strip('"')).expanduser().resolve()
    if not path.is_file():
        raise InputError('File not found: '+str(path))
    return path


def masshunter_channels(path):
    marker = False
    with path.open(encoding='utf-8-sig', newline='') as f:
        for row in csv.reader(f):
            if row and row[0].strip()=='Intensity Vs Time':
                marker = True
            if row and row[0].strip()=='Time [Sec]':
                if not marker:
                    break
                return [x.strip() for x in row[1:] if re.fullmatch(r'[A-Z][a-z]?\d+(?:\.\d+)?', x.strip())]
    raise InputError('The guided entry requires a MassHunter Intensity Vs Time export. Use a documented configuration for other formats.')


def laser_samples(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f)
        header = [x.strip() for x in next(reader)]
        if any(x not in header for x in ('Timestamp','Laser State')):
            raise InputError('The guided entry requires Timestamp and Laser State in the first header row')
        if 'Comment' not in header:
            return []
        index = header.index('Comment')
        return list(dict.fromkeys(row[index].strip() for row in reader if len(row)>index and row[index].strip()))


def _positive(value, label):
    try:
        number = float(value)
    except ValueError as exc:
        raise InputError(label+' must be a positive number') from exc
    if not math.isfinite(number) or number<=0:
        raise InputError(label+' must be a positive finite number')
    return number


def _choice(value, choices, multiple=False):
    selected = []
    for token in value.split(','):
        token = token.strip()
        if token.isdigit() and 1<=int(token)<=len(choices):
            token = choices[int(token)-1]
        if token not in choices or token in selected:
            raise InputError('Select a listed name or number, without duplicates')
        selected.append(token)
    if not multiple and len(selected)!=1:
        raise InputError('Select exactly one acquisition')
    return selected if multiple else selected[0]


def collect(ask=None, tell=None):
    ask = ask or input
    tell = tell or print
    tell(f'LA timing correction {__version__} | guided MassHunter workflow')
    tell('Choose original inputs. The tool writes a separate log only after full file re-checks.')
    ms_path = _path(ask('MassHunter signal CSV path: '))
    laser_path = _path(ask('Original laser log CSV path: '))
    channels = masshunter_channels(ms_path)
    if not channels:
        raise InputError('No supported isotope columns found')
    tell('Signal channels: '+', '.join(f'{i+1}={name}' for i,name in enumerate(channels)))
    selected = _choice(ask('Channels with useful firing contrast (names/numbers separated by commas): '), channels, True)
    names = laser_samples(laser_path)
    laser = dict(path=str(laser_path), format='laser_state_csv')
    if names:
        tell('Acquisitions: '+', '.join(f'{i+1}={name}' for i,name in enumerate(names)))
        laser['sample_selector'] = _choice(ask('Acquisition name or number: '), names)
    time_kind = ask('Laser timestamps: 1=calendar date/time, 2=numeric [1]: ').strip() or '1'
    if time_kind=='1':
        laser['time_kind'] = 'datetime'
    elif time_kind=='2':
        laser['time_kind'] = 'numeric'
        unit = ask('Laser numeric time unit (s, ms, us, min): ').strip()
        if unit not in ('s','ms','us','min'):
            raise InputError('Choose an explicit supported laser time unit')
        laser['time_unit'] = unit
    else:
        raise InputError('Choose timestamp option 1 or 2')
    kind = ask('Tolerance in 1=seconds, 2=micrometres [1]: ').strip() or '1'
    tolerance = _positive(ask('Acceptable error for this use (no universal default): '), 'Tolerance')
    analysis = {}
    if kind=='2':
        # The run converts micrometres to seconds itself, once, with the declared speed, and records the conversion.
        speed = _positive(ask('Stage speed in micrometres per second: '), 'Stage speed')
        analysis['scan_speed_um_s'] = speed
        analysis['tolerance_um'] = tolerance
        # The seconds shown here follow from the declared speed alone. The run accepts or refuses the request once it has read the log.
        stated = (f'{tolerance:g} µm ({tolerance/speed:g} s at the declared speed of {speed:g} µm s⁻¹; '
                  'the run confirms this value only after it has read the laser log)')
    elif kind!='1':
        raise InputError('Choose tolerance option 1 or 2')
    else:
        analysis['tolerance_s'] = tolerance
        stated = f'{tolerance:g} s'
    mode = ask('Export: 1=rate correction (set start offset on import), 2=rate and constant with a verified time reference [1]: ').strip() or '1'
    options = {'mode':'rate_only'}
    if mode=='2':
        tell('Supply the time of the FIRST STORED MS SAMPLE in the timestamp coordinate the receiving laser-log importer uses.')
        tell('An acquisition header time alone may omit the first sample delay. The tool cannot verify your external reference.')
        options = dict(mode='affine', ms_first_sample_at_laser_time=ask('Verified first-sample timestamp (same format/unit as the laser): ').strip(),
                       reference_note=ask('How was this time reference checked: ').strip())
    elif mode!='1':
        raise InputError('Choose export option 1 or 2')
    default = Path.cwd()/'results'/(datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')+'_'+uuid4().hex[:8])
    output = Path((ask(f'New result directory [{default}]: ').strip().strip('"') or str(default))).resolve()
    if output.exists():
        raise InputError('Result directory already exists; choose a new name')
    config = dict(ms=dict(path=str(ms_path), format='agilent_masshunter_csv', signal_columns=selected),
                  laser=laser, analysis=analysis, export=options)
    tell(f'Running diagnosis, protected export and complete CSV re-check; timing tolerance {stated}.')
    return config, str(output)
