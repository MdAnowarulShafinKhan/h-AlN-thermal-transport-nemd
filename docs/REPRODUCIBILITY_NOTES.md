# Reproducibility and curation notes

This file records important limitations of the **supplied project archive** so that public users do not over-interpret the curated repository.

## 1. Representative NEMD input is not guaranteed to be the exact production protocol

The source archive contains one representative LAMMPS script at `inputs/representative_nemd/new.langevin.lmp`. It uses a 3D periodic cell during early equilibration and contains `fix ... npt ... iso`. In LAMMPS, `iso` couples pressure control in x, y, and z. For a vacuum-separated monolayer this means the out-of-plane box can be barostatted unless the production workflow used a different input. The representative script has therefore been preserved unchanged but should be verified against the actual production inputs before it is presented as the definitive simulation protocol.

## 2. Effective thickness / cross-sectional-area convention is not encoded completely

The active numerical thickness used to convert 2D heat flow into W m^-1 K^-1 is not recoverable consistently from all supplied workbooks. The representative LAMMPS input contains only a commented `lz = 3.4` line. For porous systems, the archive also does not unambiguously identify whether gross geometric area or net solid area was used. These quantities should be made explicit in the final Methods, repository metadata, and any raw-to-kappa reproduction script.

## 3. Fifteen legacy armchair workbooks contain `_xludf` formulas

The following supplied workbooks contain Excel formula strings such as `_xludf.STDEV.S` / `_xludf.T.INV.2T`, which may display `#NAME?` in some spreadsheet engines. Their raw/hard-coded analysis values remain present, but they should ideally be regenerated with the current unified analysis script before a final archival release:

- `results/thermal_transport/armchair/pore_armchair/300K_arm_pore_12%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pore_armchair/300K_arm_pore_8%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pore_armchair/300K_arm_pore_3%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/strain_parallel_to_heat_flow_armchair/300K_arm_strain_4%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/strain_parallel_to_heat_flow_armchair/300K_arm_strain_8%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/strain_parallel_to_heat_flow_armchair/300K_arm_strain_2%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pore+strain_parallel_to_heat_flow_armchair/300K_arm_strain_8%_pore_8%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pore+strain_parallel_to_heat_flow_armchair/300K_arm_strain_4%_pore_8%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pore+strain_parallel_to_heat_flow_armchair/300K_arm_strain_0%_pore_8%_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pristine_armchair/200K_arm_pristine_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pristine_armchair/100K_arm_pristine_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pristine_armchair/300K_arm_pristine_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pristine_armchair/400K_arm_pristine_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pristine_armchair/600K_arm_pristine_full_length_analyzed.xlsx`
- `results/thermal_transport/armchair/pristine_armchair/500K_arm_pristine_full_length_analyzed.xlsx`


## 4. Nominal dopant concentrations differ among supplied structures

The supplied production-like nominal 10% B/Ga structures correspond to approximately 10.88% substitution on the cation sublattice, whereas the supplied VDOS nominal 10% structures are approximately 10.01%. If VDOS is used mechanistically to explain the transport cases, report the exact realized concentration or generate composition-matched structures.

## 5. Supplied 12% pore VDOS case is not structurally identical to the supplied 12% transport pore case

The VDOS metadata correspond to 4112 atoms versus 4576 pristine atoms (about 10.14% atom removal) and unequal Al/N counts, whereas the supplied transport-oriented 12% pore structure has exact 12% atom removal and preserved Al:N stoichiometry. Do not present the current VDOS case as a one-to-one representation of the transport structure without defining/justifying the porosity metric.

## 6. VDOS spectral resolution differs across supplied cases

Pristine, nominal 10% B/Ga and the supplied pore spectrum were processed with `block_len=4096`, while the 8% strain and 8% pore+8% strain spectra were processed with `block_len=16384`. The native frequency spacing therefore differs. For strict peak-shape comparisons, regenerate all spectra using a common processing configuration.

## 7. 10% dopant length-scaling hierarchy requires care

The supplied longest-length 10% B/Ga datasets contain multiple seeds within each realization and the seed-noise audit indicates that within-realization seed variability is important. Shorter-length workbooks do not carry the same expanded seed hierarchy. Treat infinite-length extrapolations from these mixed designs cautiously until a statistically consistent hierarchy is available at all lengths.

## 8. Raw-data completeness

The unified thermal-analysis script starts from already calculated individual thermal-conductivity values; it does not calculate kappa directly from cumulative thermostat energy, temperature gradients, geometry, and thickness. Representative raw validation data are included, but full raw energy/profile trajectories for every production condition were not present in the source archive.

## 9. VDOS raw trajectories

The processed `.dat` spectra and metadata are included, but the large `pdos.lammpstrj` velocity trajectories referenced by the metadata were not present in the source archive and therefore cannot be included here.

These notes are documentation of provenance, not a modification of the underlying numerical results.
