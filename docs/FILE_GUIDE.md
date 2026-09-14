# File and directory guide

## `results/thermal_transport/`

Contains the analyzed `.xlsx` workbooks exactly as supplied in the project archive, organized by heat-flow direction and physical perturbation. These workbooks contain individual replicate values, point statistics and, where applicable, reciprocal finite-size fits, bootstrap confidence intervals, diagnostics, sensitivity analyses, and figures embedded in the workbooks.

- `armchair/`: 36 workbooks
- `zigzag/`: 22 workbooks

## `results/temperature_bias/`

Three analyzed workbooks for the nominal thermostat driving-force / temperature-bias cases 0.02, 0.05, and 0.08 at approximately 200 nm.

## `results/vdos/`

Processed frequency-domain outputs and publication figures for pristine, 10% B, 10% Ga, 12% pore, 8% strain, and combined 8% pore + 8% strain cases. Each case may contain raw/smoothed `.dat` spectra, metadata, and PNG figures. Directional x/y/z spectra are included where they were present in the source archive.

## `validation/temperature_bias_raw/`

For each of three driving-force values and five runs (`k1`-`k5`), includes the supplied cumulative thermostat-energy file (`Ener_equ.dat`) and the project-generated energy/time and mid-region temperature-profile PNGs. The identical 4.4 MB structure file that had been copied into all 15 run folders is stored only once at `structures/pristine/AlN_width_13.7nm.lmp`.

## `validation/equilibration_sample/`

Representative temperature/energy and pressure/volume validation figure(s) and source text from the supplied 300 K, 31.6 nm example.

## `validation/steady_state_nemd/`

Representative armchair and zigzag steady-state NEMD cumulative-energy and temperature-profile data/figures.

## `structures/`

- `pristine/`: pristine/unit-cell and common thermostat-bias structure
- `doped/`: supplied equilibrated B- and Ga-substituted structures (2%, 6%, nominal 10%)
- `pores/`: supplied 3%, 8%, and 12% pore structures
- `vdos/`: supplied VDOS structures for pristine, 10% B, and 10% Ga

The source archive did not contain corresponding VDOS structure files for every VDOS condition (for example, separate 8% strain and pore+strain VDOS input structures were not available as standalone files).

## `inputs/representative_nemd/`

Contains the one representative LAMMPS input script supplied in the archive. It is archival evidence of the project workflow, not a claim that it is the exact production script for every result.

## `potentials/`

Contains the supplied `XN.tersoff` parameter file. Its comments cite the group-III-nitride Tersoff literature; see `references/REFERENCES.md` for bibliographic provenance.

## `scripts/`

- unified thermal-transport workbook analysis;
- scalar/species VDOS/PDOS analysis;
- directional x/y/z VDOS/PDOS analysis.

## `docs/` and `references/`

Curated documentation created for the GitHub package. No publisher PDFs are redistributed.
