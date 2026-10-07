"""Verify the exact files recorded in the release manifest."""
from pathlib import Path
import hashlib
import json
import sys

def main():
    root=Path(__file__).resolve().parents[1]
    try:
        manifest=json.loads((root/'FILES.json').read_text(encoding='utf-8'))
        entries=manifest['files']
        if not entries:
            raise ValueError('Manifest contains no files')
        differences=[]
        seen=set()
        for item in entries:
            rel=item['path']
            p=(root/rel).resolve()
            if p==root or root not in p.parents or rel in seen:
                raise ValueError('Invalid or repeated manifest path: '+rel)
            seen.add(rel)
            if not p.is_file():
                differences.append(dict(path=rel,reason='missing'))
            elif hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:
                differences.append(dict(path=rel,reason='sha256_mismatch'))
        print(json.dumps(dict(status='mismatch' if differences else 'verified',
            checked=len(entries),differences=differences,
            scope='Recorded bytes only; not scientific or importer validation')))
        return 2 if differences else 0
    except (OSError,ValueError,KeyError,TypeError) as exc:
        print(json.dumps(dict(status='verification_tool_error',error=str(exc))))
        return 1

if __name__=='__main__':
    raise SystemExit(main())
