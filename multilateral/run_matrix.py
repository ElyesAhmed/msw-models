#!/usr/bin/env python3
"""Run a deck with several linear solvers and print a comparison table.

Usage: ./run_matrix.py DECK.DATA [--solvers cpr_trueimpes cprw system_cpr system_cprw] [--extra "--opt=val"]
                       [--repeat 3] [--timeout 3600] [--csv results.csv] [--jobs 1]
Each solver is a --linear-solver preset name or a JSON file. Flow is started through
`mpirun -np 1` (a singleton MPI start of this build hangs). Use --jobs 1 for timings.
With --csv, one row per run is appended (see CSV_FIELDS); the table printed at the end
shows the counts and the median times over the repeats.

Every run gets its own output folder: if `out_<deck>_<config><tag>` exists already, the
run goes to `..._a2`, `..._a3` and so on, so old logs are never mixed with new ones.
The script exits with status 1 if a run failed, timed out or left no complete summary,
so that a campaign script does not take its table as done. With --allow-failed it exits
with 0 anyway (for sweeps where a failing configuration is a result); the failed runs are
still listed in the table.

--well-stats switches on the well-solve counters of the System-* JSON configurations
("well_solve_stats": "true"): the JSON file is copied into the output folder with the key
set, and flow runs with the copy. Other configurations are run unchanged.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

FLOW = os.environ.get('FLOW', 'flow')
CSV_FIELDS = ['deck', 'config', 'repeat', 'flow_version', 'flow_sha256', 'out', 'steps', 'newton',
              'linear', 'wasted_lin', 'total', 'lin_time', 'lin_setup', 'chops', 'switches',
              'well_time', 'well_inner_its', 'well_not_converged', 'well_max_rel_res',
              'well_mean_rel_res', 'rc', 'wall']
# Printed by the System solver at the end of a run with "well_solve_stats": true.
WELL_STATS = re.compile(r'System well solves \(run total\): linear_solves=\d+ applications=\d+ '
                        r'well_solves=\d+ inner_iterations=(\d+) not_converged=(\d+) time=([\d.]+) '
                        r'max_rel_residual=(\S+) mean_rel_residual=(\S+)')


def flow_sha256():
    path = shutil.which(FLOW) or FLOW
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def reserve_dir(base):
    """Create and return a new output folder: base, or base_a2, base_a3, ... if taken."""
    out, n = base, 1
    while True:
        try:
            os.makedirs(out)
            return out
        except FileExistsError:
            n += 1
            out = '%s_a%d' % (base, n)


def with_well_stats(solver, out):
    """Path of a copy of a System-* JSON config with the well-solve counters on; other solvers unchanged."""
    if not solver.endswith('.json'):
        return solver
    with open(solver) as f:
        cfg = json.load(f)
    prec = cfg.get('preconditioner', {})
    if not str(prec.get('type', '')).startswith('system'):
        return solver
    prec['well_solve_stats'] = 'true'
    path = os.path.join(out, os.path.basename(solver))
    with open(path, 'w') as f:
        json.dump(cfg, f, indent=2)
    return path


def run(deck, solver, extra, out, timeout, parallel):
    # Open MPI binds every `mpirun -np 1` to the same core, so concurrent runs
    # would share one core. Do not bind when several runs go at once; a single
    # (timing) run keeps the default binding.
    bind = ['--bind-to', 'none'] if parallel else []
    cmd = (['timeout', str(timeout)] if timeout else []) + ['mpirun', '-np', '1'] + bind + [FLOW, deck, '--output-dir=' + out, '--linear-solver=' + solver,
           '--threads-per-process=1'] + extra
    t0 = time.time()
    with open(os.path.join(out, 'run.log'), 'w') as log:
        rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='1'))
    return rc, time.time() - t0


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
    m = WELL_STATS.search(txt)
    if m:
        r.update(well_inner_its=int(m.group(1)), well_not_converged=int(m.group(2)),
                 well_time=float(m.group(3)), well_max_rel_res=float(m.group(4)),
                 well_mean_rel_res=float(m.group(5)))
    return r


def complete(r):
    return all(not (isinstance(r[k], float) and math.isnan(r[k])) for k in ('steps', 'newton', 'linear', 'total'))


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
    p.add_argument('--well-stats', action='store_true',
                   help='switch on the well-solve counters of System-* JSON configs')
    p.add_argument('--allow-failed', action='store_true',
                   help='exit with 0 even if runs failed (they are still listed)')
    a = p.parse_args()
    sha = flow_sha256()
    name = lambda s: os.path.splitext(os.path.basename(s))[0]
    deckname = os.path.splitext(os.path.basename(a.deck))[0]
    tag = lambda r: a.tag + ('_r%d' % r if a.repeat > 1 else '')
    # Reserve the output folders up front, in order, so parallel runs cannot collide.
    jobs = [(s, rep, reserve_dir('out_%s_%s%s' % (deckname, name(s), tag(rep))))
            for s in a.solvers for rep in range(a.repeat)]
    used = [with_well_stats(s, out) if a.well_stats else s for s, rep, out in jobs]
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        res = list(ex.map(lambda j: run(a.deck, j[1], a.extra.split(), j[0][2], a.timeout, a.jobs > 1),
                          zip(jobs, used)))
    rows = []
    for (solver, rep, out), (rc, wall) in zip(jobs, res):
        r = parse(out)
        r.update(deck=os.path.basename(a.deck), config=name(solver), repeat=rep, rc=rc,
                 wall=round(wall, 2), flow_sha256=sha, out=out)
        rows.append(r)
    if a.csv:
        fields = CSV_FIELDS
        if os.path.exists(a.csv):
            with open(a.csv, newline='') as f:
                fields = next(csv.reader(f), CSV_FIELDS)
        with open(a.csv, 'a', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
            if not os.path.getsize(a.csv):
                w.writeheader()
            w.writerows(rows)
    print('%-22s %5s %6s %6s %6s %8s %8s %8s %6s %6s %3s' % (
        'solver', 'steps', 'newton', 'linear', 'wst%', 'total s', 'lin s', 'setup s', 'chops', 'switch', 'rc'))
    problems = []
    for solver in a.solvers:
        rs = [r for r in rows if r['config'] == name(solver)]
        ok = [r for r in rs if r['rc'] == 0 and complete(r)]
        for r in rs:
            if r not in ok:
                problems.append('%s repeat %d: rc=%d%s (%s)' % (
                    r['config'], r['repeat'], r['rc'], '' if complete(r) else ', no complete summary', r['out']))
        if not ok:
            print('%-22s %s' % (name(solver), 'no successful run'))
            continue
        counts = {(r['steps'], r['newton'], r['linear']) for r in ok}
        if len(counts) > 1:
            problems.append('%s: counts differ between repeats: %s' % (name(solver), sorted(counts)))
        med = lambda k: statistics.median(r[k] for r in ok)
        print('%-22s %5s %6s %6s %6.1f %8.2f %8.2f %8.2f %6d %6d %3d' % (
            ok[0]['config'], ok[0]['steps'], ok[0]['newton'], ok[0]['linear'], ok[0]['wasted_lin'],
            med('total'), med('lin_time'), med('lin_setup'), ok[0]['chops'], ok[0]['switches'],
            max(r['rc'] for r in rs)))
        if 'well_time' in ok[0]:
            print('%-22s well solves: %.2f s, %d inner its, %d not converged, max rel res %.2e, mean %.2e' % (
                '', med('well_time'), ok[0]['well_inner_its'], ok[0]['well_not_converged'],
                ok[0]['well_max_rel_res'], ok[0]['well_mean_rel_res']))
    print('flow %s sha256 %s' % (shutil.which(FLOW) or FLOW, sha))
    for msg in problems:
        print('WARNING: ' + msg)
    failed = any(r['rc'] != 0 or not complete(r) for r in rows)
    sys.exit(1 if failed and not a.allow_failed else 0)


if __name__ == '__main__':
    main()
