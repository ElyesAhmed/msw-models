#!/usr/bin/env python3
"""Generate a multilateral multisegment-well stress case on the sector01 models.

The reservoir description (grid, properties, PVT, saturation functions, equilibration) and
the water injector are taken unchanged from the base deck of the chosen model:
  sector01     ../sector01/SECTOR01_MSW1.DATA                 (15 x 25 x 134, injector A3H)
  sector01x20  ../sector01x20/SECTOR01X20_MSW2.DATA          (15 x 500 x 134, injector A3HX20,
               the sector01 grid refined 20 times along J; needs the Git LFS data)
The single-branch producers of the base deck are replaced by `--wells` fishbone producers.
Each producer has a horizontal mother bore along J at a fixed (I, K) in the oil column and
`--laterals` short laterals along I. The well layout is defined on the sector01 J index and
scaled to the refined grid, and the mother bore keeps the same number of segments on both
grids, so the case stresses the number of wells and branches rather than segments.

Usage:
  ./make_multilateral.py [--model sector01|sector01x20] [--model-dir DIR] [--wells 6]
                         [--laterals 3] [--cells-per-segment 2] [--months 12]
                         [--orat 1200] [--inj-rate 8000] [--out SECTOR01_MULTILAT.DATA]
"""
import argparse
import os
import re

import numpy as np

import grdecl

HERE = os.path.dirname(os.path.abspath(__file__))

MODELS = {
    'sector01': dict(
        dir=os.path.join(HERE, '..', 'sector01'),
        deck='SECTOR01_MSW1.DATA', actnum='ACTNUM_BLOCK.GRDECL',
        dims=(15, 25, 134), jscale=1, injector='A3H'),
    'sector01x20': dict(
        dir=os.path.join(HERE, '..', 'sector01x20'),
        deck='SECTOR01X20_MSW2.DATA', actnum='GRID.GRDECL',  # ACTNUM.GRDECL is empty here
        dims=(15, 500, 134), jscale=20, injector='A3HX20'),
}
GOC, OWC = 1542.5, 1568.5
DX, DY = 150.0, 150.0  # approximate sector01 cell sizes (m) along I and J

# Producer slots on the sector01 J index: (I, K, J of the laterals, lateral direction along I).
# Columns I=5 and I=11 at three depths, then I=8 at two; laterals are staggered in J so that
# laterals of neighbouring wells never share a cell.
SLOTS = [
    (5, 62, (4, 10, 16), (+1, -1, +1)),
    (11, 62, (7, 13, 19), (-1, +1, -1)),
    (5, 66, (4, 10, 16), (+1, -1, +1)),
    (11, 66, (7, 13, 19), (-1, +1, -1)),
    (5, 70, (4, 10, 16), (+1, -1, +1)),
    (11, 70, (7, 13, 19), (-1, +1, -1)),
    (8, 64, (5, 12, 18), (+1, -1, +1)),
    (8, 68, (5, 12, 18), (-1, +1, -1)),
]
MOTHER_J = (2, 20)
LATERAL_LEN = 3  # cells along I


def load_props(model):
    m = MODELS[model]
    nx, ny, nz = m['dims']
    inc = os.path.join(m['dir'], 'INCLUDE') + '/'
    for f in ('GRID.GRDECL', m['actnum']):
        with open(inc + f) as fh:
            if fh.readline().startswith('version https://git-lfs'):
                raise SystemExit('%s is a Git LFS pointer; run "git lfs pull" in the model repository' % (inc + f))
    act = grdecl.read_keyword(inc + m['actnum'], 'ACTNUM').reshape(nz, ny, nx)
    depth = grdecl.cell_depths(inc + 'GRID.GRDECL', nx, ny, nz)
    return act, depth


def fine_j(j, s):
    """J index range on the refined grid covered by sector01 cell j."""
    return range((j - 1) * s + 1, j * s + 1)


def build_well(name, slot, nlat, cps, model, act, depth, used):
    """Return the deck text and a geometry record for one fishbone producer."""
    nx, ny, _ = MODELS[model]['dims']
    s = MODELS[model]['jscale']
    i0, k0, lat_j, lat_dir = slot
    conns = []      # (i, j, k, dir, branch, seg)
    segs = []       # (seg, branch, outlet, md, tvd)
    geom = {'name': name, 'mother': [], 'laterals': []}

    def ok(i, j, k):
        return 1 <= i <= nx and 1 <= j <= ny and act[k - 1, j - 1, i - 1] > 0 and (i, j, k) not in used

    # Mother bore (branch 1): segment 1 is the top segment at the heel; one segment per
    # `cps` sector01 cells, i.e. per cps * jscale cells of the refined grid.
    j_heel = fine_j(MOTHER_J[0], s)[0]
    top_tvd = depth[k0 - 1, j_heel - 1, i0 - 1]
    md, seg = 100.0, 1
    junction = {}
    for n, jc in enumerate(range(MOTHER_J[0], MOTHER_J[1] + 1)):
        if n % cps == 0:
            seg += 1
            md += cps * DY
            segs.append((seg, 1, seg - 1, md, depth[k0 - 1, fine_j(jc, s)[0] - 1, i0 - 1]))
        junction[jc] = seg
        for j in fine_j(jc, s):
            geom['mother'].append((i0, j, k0))
            if ok(i0, j, k0):
                used.add((i0, j, k0))
                conns.append((i0, j, k0, 'Y', 1, seg))
    # Laterals (branches 2..) leave the mother bore along I, from the middle of the sector01 cell.
    branch = 1
    for jc, d in list(zip(lat_j, lat_dir))[:nlat]:
        branch += 1
        outlet = junction[jc]
        lmd = [x for x in segs if x[0] == outlet][0][3]
        j = fine_j(jc, s)[s // 2]
        lat = []
        for n in range(LATERAL_LEN):
            i = i0 + d * (n + 1)
            if not 1 <= i <= nx:
                break
            if n % cps == 0:
                seg += 1
                lmd += cps * DX
                segs.append((seg, branch, outlet, lmd, depth[k0 - 1, j - 1, i - 1]))
                outlet = seg
            lat.append((i, j, k0))
            if ok(i, j, k0):
                used.add((i, j, k0))
                conns.append((i, j, k0, 'X', branch, seg))
        geom['laterals'].append(lat)

    w = []
    w.append('WELSPECS\n%-6s PROD %3d %3d %.1f OIL 1* 1* SHUT 1* 1* 1* 1* /\n/\n'
             % (name, i0, j_heel, top_tvd))
    w.append('COMPORD\n%-6s INPUT /\n/\n' % name)
    w.append('COMPDAT\n' + ''.join('%-6s %3d %3d %3d %3d OPEN 1* 1* 0.2159 1* 0 1* %s /\n'
                                   % (name, i, j, k, k, dr) for i, j, k, dr, _, _ in conns) + '/\n')
    w.append('WELSEGS\n%-6s %.2f 100.0 1* ABS HF- HO /\n' % (name, top_tvd)
             + ''.join('%4d %4d %3d %4d %9.2f %8.2f 0.159 0.00065 /\n' % (sg, sg, b, o, m, t)
                       for sg, b, o, m, t in segs) + '/\n')
    # COMPSEGS: each connection gets its segment explicitly, through a measured-depth interval
    # inside that segment.
    segmd = {sg: m for sg, _, _, m, _ in segs}
    w.append('COMPSEGS\n%-6s /\n' % name
             + ''.join('%3d %3d %3d %3d %9.2f %9.2f 1* 1* 1* 1* %d /\n'
                       % (i, j, k, b, segmd[sg] - 1.0, segmd[sg], sg) for i, j, k, _, b, sg in conns)
             + '/\n')
    geom.update(nconn=len(conns), nseg=len(segs), nbranch=branch)
    return ''.join(w), geom


def base_sections(model, model_dir):
    """RUNSPEC to SUMMARY of the base deck, and the injector definition."""
    m = MODELS[model]
    txt = open(os.path.join(model_dir, m['deck'])).read()
    head = txt[:txt.index('\nSCHEDULE')]
    sched = txt[txt.index('\nSCHEDULE'):]
    inc = re.search(r"INCLUDE\s*\n\s*'(INCLUDE/WELL_%s\.INC)'\s*/" % m['injector'], sched)
    if inc:
        inj = "INCLUDE\n  '%s/%s' /\n" % (os.path.relpath(model_dir), inc.group(1))
    else:
        inj = sched[sched.index('WELSPECS\n%s ' % m['injector']):sched.index('WCONPROD')]
    return head, inj


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', choices=sorted(MODELS), default='sector01')
    p.add_argument('--model-dir', help='directory of the base deck (default: see MODELS)')
    p.add_argument('--wells', type=int, default=6)
    p.add_argument('--laterals', type=int, default=3)
    p.add_argument('--cells-per-segment', type=int, default=2,
                   help='sector01 cells per segment (the mother bore keeps this on the refined grid)')
    p.add_argument('--months', type=int, default=12)
    p.add_argument('--orat', type=float, default=1200.0)
    p.add_argument('--inj-rate', type=float, default=8000.0)
    p.add_argument('--out', default='SECTOR01_MULTILAT.DATA')
    a = p.parse_args()
    assert 1 <= a.wells <= len(SLOTS)
    assert 1 <= a.laterals <= 3
    model_dir = a.model_dir or MODELS[a.model]['dir']
    MODELS[a.model]['dir'] = model_dir

    act, depth = load_props(a.model)
    head, inj = base_sections(a.model, model_dir)
    s = MODELS[a.model]['jscale']
    max_conn = (MOTHER_J[1] - MOTHER_J[0] + 1) * s + a.laterals * LATERAL_LEN
    head = re.sub(r'WELLDIMS\s*\n[^/]*/', 'WELLDIMS\n  %d %d 2 %d /' % (a.wells + 1, max(400, max_conn), a.wells + 1), head)
    head = re.sub(r'WSEGDIMS\s*\n[^/]*/', 'WSEGDIMS\n  %d 60 %d /' % (a.wells + 1, a.laterals + 2), head)
    # Relative to the current directory, where the deck is written.
    head = head.replace("'INCLUDE/", "'%s/INCLUDE/" % os.path.relpath(model_dir))

    used, wells, geoms = set(), [], []
    for n in range(a.wells):
        txt, g = build_well('P%02d' % (n + 1), SLOTS[n], a.laterals, a.cells_per_segment,
                            a.model, act, depth, used)
        wells.append(txt)
        geoms.append(g)

    names = [g['name'] for g in geoms]
    sched = ['\nSCHEDULE\nRPTRST\n BASIC=3 FREQ=3 NORST=1 /\n'] + wells + [inj]
    sched.append('WCONPROD\n' + ''.join('%-6s OPEN ORAT %.0f 1* 1* 1* 1* 100.0 1* 1* 1* /\n' % (w, a.orat)
                                        for w in names) + '/\n')
    sched.append('WCONINJE\n%s WAT OPEN RATE %.0f 1* 200. 1* 1* /\n/\n' % (MODELS[a.model]['injector'], a.inj_rate))
    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    dates = ['  1 %s %d /\n' % (months[(m + 1) % 12], 2000 + (m + 1) // 12) for m in range(a.months)]
    sched.append('DATES\n' + ''.join(dates) + '/\nEND\n')

    with open(a.out, 'w') as f:
        f.write('-- Generated by make_multilateral.py from %s (ODbL).\n' % MODELS[a.model]['deck'])
        f.write('-- %s\n' % vars(a))
        f.write(head + ''.join(sched))
    np.save(os.path.splitext(a.out)[0] + '_geom.npy', np.array(geoms, dtype=object), allow_pickle=True)
    tot = lambda key: sum(g[key] for g in geoms)
    print('%s: %d producers, %d branches, %d segments, %d connections'
          % (a.out, len(geoms), tot('nbranch'), tot('nseg'), tot('nconn')))


if __name__ == '__main__':
    main()
