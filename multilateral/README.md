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
- `cpr`, `cprw`, `system_cpr` and `system_cprw` are the presets written out. `cprw` is OPM's default (`--linear-solver=cprw`, BiCGStab outer solver); include it in every comparison.
- `system_cprw_equiv` is System-CPRW with the coarse diagonal of classic CPRW (`well_coarse_diagonal: auto`). With exact well solves it gives the same iterates as `cprw`.
- `*_noid` keep the contracted coarse equation of pressure-controlled wells instead of a trivial one (`well_identity_on_pressure_control: false`).
- `cprw_contract_d*` contract the well block D for the coarse diagonal of multisegment wells in classic CPRW (`well_coarse_diagonal: contract_d`), as `system_cprw` does.
- `*_fgmres` use an FGMRES outer solver.
- `*_inexact_*` use iterative inner well solves (FGMRES + DILU). The `*_inexact_it5_{classic,no_prolongation,full}` files differ from `*_inexact_it5` only in `well_transfer`, and `system_cprw_inexact_it1_{no_prolongation,full}` differ from `system_cprw_inexact_it1` only in `well_transfer`.

Classic and System configurations with the same coarse matrix, which give the same iterates with exact well solves:

| Classic CPRW | System-CPRW | Coarse diagonal (MSW) | Pressure-controlled wells |
|---|---|---|---|
| `cprw` | `system_cprw_equiv` | row sum | trivial equation |
| `cprw_noid` | `system_cprw_equiv_noid` | row sum | contracted equation |
| `cprw_contract_d` | `system_cprw` | contract D | trivial equation |
| `cprw_contract_d_noid` | `system_cprw_noid` | contract D | contracted equation |

Pass the JSON files rather than the presets, and always set `well_identity_on_pressure_control`: a System JSON file without it gets `false`, a classic CPRW one gets `true`, and `--linear-solver=system_cprw` sets `true`.

## OPM build for the paper

Build flow from branch `cg-paper-system-cprw` of https://github.com/ElyesAhmed/opm-simulators, with opm-common `7d592ae23` and opm-grid `3a6406b34` (both on OPM master of 10 Nov 2026; the same commits are also on branch `cg-paper-system-cprw` of https://github.com/ElyesAhmed/opm-common and https://github.com/ElyesAhmed/opm-grid). The branch is OPM master `04f48b499` plus:

1. **AMG reuse fix:** the System solver follows `--cpr-reuse-setup`, as classic CPRW does. Without it, System-* keep the AMG hierarchy of the first linear system for the whole run, and their results are not valid.
2. **`well_identity_on_pressure_control` for classic CPRW,** needed by `cprw_noid` and `cprw_contract_d_noid`. An older build ignores the key without a warning and runs plain `cprw`, so check the build as below.
3. **`well_solve_stats`** (System-* only, `"preconditioner": {"well_solve_stats": "true"}`): at the end of the run flow prints `System well solves (run total): ...` with the time spent in well solves, the inner iterations, the number of inner solves that did not converge and the true relative residual of the inner solves. `run_matrix.py` reads it. Leave it off for timing runs.

### Check the build

Counts only, so the runs may go in parallel (`--jobs`). With the build above, the counts must be exactly these (time steps / Newton / linear iterations):

| Deck | Configurations | Expected |
|---|---|---|
| `opm-tests/spe9/SPE9.DATA` | `cprw`, `cprw_contract_d`, `system_cprw_equiv`, `system_cprw` | 93 / 237 / 392 |
| `opm-tests/spe9/SPE9.DATA` | `cprw_noid`, `cprw_contract_d_noid`, `system_cprw_equiv_noid`, `system_cprw_noid` | 93 / 237 / 427 |
| `../sector01x20/SECTOR01X20_MSW2.DATA` | `cprw`, `system_cprw_equiv` | 42 / 307 / 783 |
| `../sector01x20/SECTOR01X20_MSW2.DATA` | `cprw_contract_d`, `system_cprw` | 42 / 307 / 782 |
| `SECTOR01_ML8_S2.DATA` | `cprw` | 99 / 818 / 1682 |

```bash
FLOW=/path/to/flow ./run_matrix.py /path/to/opm-tests/spe9/SPE9.DATA --jobs 8 \
  --solvers configs/cprw.json configs/cprw_noid.json configs/cprw_contract_d.json configs/cprw_contract_d_noid.json \
            configs/system_cprw_equiv.json configs/system_cprw_equiv_noid.json configs/system_cprw.json configs/system_cprw_noid.json
```

If `cprw_noid` gives 392 instead of 427, the build does not have the classic CPRW option. If `system_cprw_equiv` gives 45 / 320 / 1342 on Sector-MSW, the build does not have the AMG reuse fix.

On the fishbone decks the matched classic and System pairs agree only to a few per cent, because the many time-step cuts make the Newton path sensitive to rounding. On `SECTOR01_ML8_S2` the four `*_noid` configurations all give 899–919 linear iterations against 1682–1784 with the trivial equation: keeping the contracted equation of pressure-controlled wells helps classic CPRW as much as System-CPRW.

`run_matrix.py` starts flow with `mpirun -np 1`. When it runs several at once (`--jobs` > 1) it adds `--bind-to none`, since Open MPI otherwise binds every run to the same core. It exits with status 1 if a run failed, timed out or left no complete summary, and gives every run a new output folder.

`run_matrix.py DECK --solvers cfg.json ... [--repeat 3] [--timeout 3600] [--csv out.csv] [--jobs N] [--well-stats] [--allow-failed]` runs several configurations and prints iteration and timing medians. `FLOW=/path/to/flow` selects the binary. Run timing runs one at a time.
- `--well-stats` turns on the well-solve counters of the System-* JSON files: each file is copied into the run's output folder with `"well_solve_stats": "true"`. The other configurations run unchanged.
- `--allow-failed` exits with status 0 even if runs failed, for sweeps where a failing configuration is a result. The failed runs are still listed.

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
