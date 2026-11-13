"""Minimal GRDECL reader for the sector01 corner-point grid (no external dependencies besides numpy)."""
import re
import numpy as np


def read_keyword(path, kw):
    """Return the values of keyword `kw` in a GRDECL file as a float array (expands n*v)."""
    vals, on = [], False
    with open(path) as f:
        for line in f:
            line = line.split('--')[0].strip()
            if not line:
                continue
            if not on:
                if re.fullmatch(kw, line.split()[0]):
                    on = True
                continue
            end = '/' in line
            for tok in line.split('/')[0].split():
                if '*' in tok:
                    n, v = tok.split('*')
                    vals.extend([float(v)] * int(n))
                else:
                    vals.append(float(tok))
            if end:
                break
    return np.array(vals)


def cell_depths(grid_path, nx, ny, nz):
    """Cell-centre depths (nz, ny, nx) from ZCORN: mean of the 8 corner depths."""
    z = read_keyword(grid_path, 'ZCORN').reshape(2 * nz, 2 * ny, 2 * nx)
    z = z.reshape(nz, 2, ny, 2, nx, 2).mean(axis=(1, 3, 5))
    return z
