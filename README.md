# Thermal transport in monolayer h-AlN: molecular-dynamics data and analysis

This repository is a curated research package for classical molecular-dynamics studies of heat transport in monolayer hexagonal aluminum nitride (h-AlN). It contains the project-owned analysis scripts, analyzed Excel workbooks, representative LAMMPS input/data files, interatomic-potential file, VDOS/PDOS post-processing outputs, figures, and representative raw NEMD validation data supplied in the source project archive.

## Study dimensions represented here

- pristine armchair and zigzag thermal transport from 100-600 K;
- finite-length extrapolation and width checks;
- uniaxial-strain cases;
- pore/porosity cases;
- combined pore + strain cases;
- B- and Ga-substitution cases and selected dopant + strain screening;
- thermostat driving-force / temperature-bias sensitivity checks;
- VDOS/PDOS and directional VDOS analysis;
- representative equilibration, cumulative thermostat-energy, and temperature-profile validation outputs.

## Repository layout

```text
scripts/                    Python analysis and VDOS/PDOS scripts
inputs/                     Representative LAMMPS input script
potentials/                 Tersoff parameter file
structures/                 Pristine, doped, porous, and VDOS LAMMPS structures
results/thermal_transport/  58 analyzed transport workbooks
results/temperature_bias/   3 analyzed thermostat-bias workbooks
results/vdos/               Processed VDOS/PDOS data, metadata, and figures
validation/                 Representative raw/figure validation outputs
docs/                       File guide, reproducibility notes, exclusions, benchmark table
references/                 Bibliographic links; publisher PDFs are not redistributed
MANIFEST.csv                File-level provenance, sizes, and SHA-256 hashes
SHA256SUMS.txt               Integrity hashes for the curated repository files
```

The repository contains **all 61 analyzed Excel workbooks** from the project archive. Fifteen identical copies of the same 13.7 nm-width LAMMPS structure used in the thermostat-bias runs were deduplicated into `structures/pristine/AlN_width_13.7nm.lmp`.

## Python environment

Python 3.10+ is recommended.

```bash
python -m venv .venv
source .venv/bin/activate        # Linux/macOS/WSL
# .venv\Scripts\activate       # Windows PowerShell/cmd variant
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### Thermal-transport workbook analysis

```bash
python scripts/thermal_transport/thermal_transport_analysis_ALL_CASES_FINAL.py input.xlsx
```

See `scripts/thermal_transport/ORIGINAL_SCRIPT_README.txt` for supported workbook layouts, bootstrap settings, finite-size extrapolation rules, and optional arguments.

### VDOS/PDOS analysis

The VDOS scripts require a LAMMPS trajectory containing at least `id type vx vy vz`.

Example for pristine h-AlN using the settings stored in the supplied metadata:

```bash
python scripts/vdos/pdos_publication.py \
  --dump pdos.lammpstrj \
  --timestep_fs 0.5 \
  --dump_interval_steps 1 \
  --stride 2 \
  --block_len 4096 \
  --species_map "1:Al,2:N" \
  --mass_map "Al=26.9815385,N=14.0067" \
  --primitive_atoms 2 \
  --smooth_sigma_thz 0.38 \
  --xmin 0 --xmax 45 \
  --prefix PDOS_hAlN
```

For Cartesian x/y/z-resolved spectra, use `scripts/vdos/pdos_publication_directional.py` with the same trajectory metadata and a structure whose Cartesian axes correspond to the intended crystallographic directions.

## LAMMPS inputs

`inputs/representative_nemd/new.langevin.lmp` is the representative input script present in the supplied archive; it is preserved rather than silently rewritten. `structures/` and `potentials/XN.tersoff` contain the supplied data/potential files.

**Important:** read `docs/REPRODUCIBILITY_NOTES.md` before treating the representative input as the exact production protocol. The source archive does not contain exact production inputs/raw trajectories for every reported condition.

## Data integrity and provenance

`MANIFEST.csv` records the source-archive path, curated repository path, file size, SHA-256 digest, and file role. `SHA256SUMS.txt` provides repository integrity hashes.

Third-party literature PDFs, journal author guidelines, and the prior ChatGPT transcript were intentionally excluded. Their relevant bibliographic information is summarized in `references/REFERENCES.md` and the exclusions are documented in `docs/EXCLUDED_FILES.md`.

## Reproducibility status

This repository preserves the project archive faithfully but is **not yet an end-to-end raw-trajectory-to-final-paper reproduction package**. The main limitations and known metadata/statistical issues identified during curation are listed transparently in `docs/REPRODUCIBILITY_NOTES.md`.

## Citation

GitHub can read the included `CITATION.cff`. Update the author list and replace the repository citation with the final article citation/DOI when the associated manuscript is published.

## License

No open-source or open-data license has yet been assigned. See LICENSE_NOTICE.md for the current licensing status.
