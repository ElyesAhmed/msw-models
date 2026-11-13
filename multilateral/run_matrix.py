#!/usr/bin/env python3
"""Run a deck with several linear solvers and print a comparison table.

Usage: ./run_matrix.py DECK.DATA [--solvers cpr_trueimpes cprw system_cpr system_cprw] [--extra "--opt=val"]
                       [--repeat 3] [--timeout 3600] [--csv results.csv] [--jobs 1]
Each solver is a --linear-solver preset name or a JSON file. Flow is started through
`mpirun -np 1` (a singleton MPI start of this build hangs). Use --jobs 1 for timings.
With --csv, one row per run is appended (see CSV_FIELDS); the table printed at the end
shows the median over the repeats.
"""
import csv
import statistics
import argparse
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

FLOW = os.environ.get('FLOW', 'flow')
CSV_FIELDS = ['deck', 'config', 'repeat', 'flow_version', 'steps', 'newton', 'linear', 'wasted_lin',
              'total', 'lin_time', 'lin_setup', 'chops', 'switches', 'rc', 'wall']


def run(deck, solver, extra, tag, timeout):
    out = 'out_%s_%s%s' % (os.path.splitext(os.path.basename(deck))[0],
                           os.path.splitext(os.path.basename(solver))[0], tag)
    os.makedirs(out, exist_ok=True)
    cmd = (['timeout', str(timeout)] if timeout else []) + ['mpirun', '-np', '1', FLOW, deck, '--output-dir=' + out, '--linear-solver=' + solver,
           '--threads-per-process=1'] + extra
    t0 = time.time()
    with open(os.path.join(out, 'run.log'), 'w') as log:
        rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='1'))
    return solver, out, rc, time.time() - t0


def parse(out):
    txt = open(os.path.join(out, 'run.log')).read()
    prt = ''.join(open(os.path.join(out, f)).read() for f in os.listdir(out) if f.endswith('.PRT'))
    def grab(pat, cast=float):
        m = re.search(pat, txt)
        return cast(m.group(1)) if m else float('nan')
    r = {
        'steps': grab(r'Number of timesteps:\s+(\d+)', int),
        'newton': grab(r'Overall Newton Iterations:\s+(\d+)', int),
        'linear': grab(r'Overall Linear Iterations:\s+(\d+)', int),
        'wasted_lin': grab(r'Overall Linear Iterations:.*?;\s*([\d.]+)%'),
        'total': grab(r'Simulation time:\s+([\d.]+)'),
        'lin_time': grab(r'Linear solve time:\s+([\d.]+)'),
        'lin_setup': grab(r'Linear setup:\s+([\d.]+)'),
        'chops': len(re.findall(r'(?i)chopp', prt)),
        'switches': len(re.findall(r'control mode changed', prt)),
        'flow_version': (re.search(r'(flow \S+|OPM Flow version \S+)', txt) or [''])[0],
    }
    return r


def main():
    p = argparse.ArgumentParser()
    p.add_argument('deck')
    p.add_argument('--solvers', nargs='+', default=['cpr_trueimpes', 'cprw', 'system_cpr', 'system_cprw'])
    p.add_argument('--extra', default='')
    p.add_argument('--tag', default='')
    p.add_argument('--jobs', type=int, default=1, help='concurrent runs (use 1 for timings)')
    p.add_argument('--repeat', type=int, default=1)
    p.add_argument('--timeout', type=int, default=0, help='seconds per run, 0 for none')
    p.add_argument('--csv', help='append one row per run to this CSV file')
    a = p.parse_args()
    jobs = [(s, r) for s in a.solvers for r in range(a.repeat)]
    tag = lambda r: a.tag + ('_r%d' % r if a.repeat > 1 else '')
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        res = list(ex.map(lambda j: run(a.deck, j[0], a.extra.split(), tag(j[1]), a.timeout) + (j[1],), jobs))
    rows = []
    for solver, out, rc, wall, rep in res:
        r = parse(out)
        r.update(deck=os.path.basename(a.deck), config=os.path.splitext(os.path.basename(solver))[0],
                 repeat=rep, rc=rc, wall=round(wall, 2))
        rows.append(r)
    if a.csv:
        new = not os.path.exists(a.csv)
        with open(a.csv, 'a', newline='') as f:
            w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction='ignore')
            if new:
                w.writeheader()
            w.writerows(rows)
    print('%-22s %5s %6s %6s %6s %8s %8s %8s %6s %6s %3s' % (
        'solver', 'steps', 'newton', 'linear', 'wst%', 'total s', 'lin s', 'setup s', 'chops', 'switch', 'rc'))
    for solver in a.solvers:
        rs = [r for r in rows if r['config'] == os.path.splitext(os.path.basename(solver))[0]]
        med = lambda k: statistics.median(r[k] for r in rs)
        print('%-22s %5s %6s %6s %6.1f %8.2f %8.2f %8.2f %6d %6d %3d' % (
            rs[0]['config'], rs[0]['steps'], rs[0]['newton'], rs[0]['linear'], rs[0]['wasted_lin'],
            med('total'), med('lin_time'), med('lin_setup'), rs[0]['chops'], rs[0]['switches'],
            max(r['rc'] for r in rs)))


if __name__ == '__main__':
    main()
