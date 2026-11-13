# Fishbone multilateral stress cases

These decks are derived from `../sector01` and `../sector01x20`. The reservoir description (grid, properties, PVT, saturation functions, equilibration) and the water injector are unchanged; the single-branch producers are replaced by fishbone multisegment producers. Each producer has a horizontal mother bore along J and three laterals along I. The producers are on oil-rate control (3000 Sm3/day) with a 100 bar BHP limit and switch to BHP control during the run.

They are used in the System-CPRW paper (Computational Geosciences, ECMOR 2026 collection) to test CPRW-type preconditioners on many coupled multisegment wells.

| Deck | Grid | Wells / branches / segments / connections | Period |
|---|---|---|---|
| `SECTOR01_ML8_S2` | 15 x 25 x 134 | 8 / 32 / 128 / 209 | 24 months |
| `SECTOR01_ML8_S1` | 15 x 25 x 134 | 8 / 32 / 224 / 209 | 24 months |
| `SECTOR01X20_ML8_S2` | 15 x 500 x 134 (about 280k cells) | 8 / 32 / 128 / 2964 | 24 months |
| `SECTOR01X20_ML8_S1` | 15 x 500 x 134 | 8 / 32 / 224 / 2964 | 24 months |
| `SECTOR01X20_ML8_S2_6M` | 15 x 500 x 134 | as `_S2` | 6 months |
| `SECTOR01X20_ML8_BHP_1M` | 15 x 500 x 134 | as `_S2`, all producers on BHP control from the start | 1 month |

S2 means two cells per segment; S1 means one. On the x20 grid the mother bore keeps the same number of segments, so the refined decks stress the number of connections rather than segments.

The decks include their grid and property files from `../sector01/INCLUDE` and `../sector01x20/INCLUDE` by relative path, so keep this folder inside the repository. The x20 files are Git LFS objects; run `git lfs pull` first.

## Running

```bash
mpirun -np 1 flow SECTOR01_ML8_S2.DATA --linear-solver=configs/cprw.json --threads-per-process=1
```

`configs/` holds the linear-solver JSON files of the paper:
- `cpr`, `cprw`, `system_cpr` and `system_cprw` are the presets written out.
- `system_cprw_noid` has no trivial coarse equation for pressure-controlled wells.
- `*_fgmres` use an FGMRES outer solver.
- `*_inexact_*` use iterative inner well solves.

Pass the JSON files rather than the presets. A JSON file without `well_identity_on_pressure_control` gets `false`, while `--linear-solver=system_cprw` sets it to `true`.

The System-* configurations need an opm-simulators build in which the System solver follows `--cpr-reuse-setup`, as classic CPRW does. Without that fix, the System solver keeps the AMG hierarchy of the first linear system for the whole run.

`run_matrix.py DECK --solvers cfg.json ... [--repeat 3] [--timeout 3600] [--csv out.csv] [--jobs N]` runs several configurations and prints iteration and timing medians. `FLOW=/path/to/flow` selects the binary. Run timing runs one at a time.

## Regenerating

Run these from this folder; `plot_model.py DECK` plots the wells.

```bash
./make_multilateral.py --wells 8 --orat 3000 --months 24 --out SECTOR01_ML8_S2.DATA
./make_multilateral.py --wells 8 --orat 3000 --months 24 --cells-per-segment 1 --out SECTOR01_ML8_S1.DATA
./make_multilateral.py --model sector01x20 --wells 8 --orat 3000 --months 24 --out SECTOR01X20_ML8_S2.DATA
./make_multilateral.py --model sector01x20 --wells 8 --orat 3000 --months 24 --cells-per-segment 1 --out SECTOR01X20_ML8_S1.DATA
./make_multilateral.py --model sector01x20 --wells 8 --orat 3000 --months 6 --out SECTOR01X20_ML8_S2_6M.DATA
```

`SECTOR01X20_ML8_BHP_1M` is `_S2_6M` edited by hand: `WCONPROD` set to BHP 100 bar for all producers, and the schedule cut to the first month.

The decks are derived from Equinor data and keep its licence (Open Database License, see the deck headers).
