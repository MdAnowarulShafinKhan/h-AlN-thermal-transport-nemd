# Scripts

## Thermal transport

`thermal_transport/thermal_transport_analysis_ALL_CASES_FINAL.py` analyzes already-calculated individual thermal conductivities in Excel input layouts. It can perform point statistics, Student-t confidence intervals, reciprocal finite-size fits, bootstrap uncertainty, sensitivity diagnostics, hierarchical realization/seed treatment, and special checks for width/delta/direction series.

It **does not** calculate kappa directly from raw NEMD heat-flux and temperature-gradient files.

## VDOS/PDOS

`vdos/pdos_publication.py` computes MD-derived velocity power spectra / VDOS/PDOS from a LAMMPS velocity dump.

`vdos/pdos_publication_directional.py` additionally resolves Cartesian x/y/z contributions while enforcing total = x + y + z with one common normalization.

The scripts are supplied as research code and should be used with the metadata/provenance notes in this repository.
