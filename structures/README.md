# Structures

This directory collects the standalone LAMMPS data/structure files available in the source project archive.

- `pristine/AlN_paper_like_unit_cell.lmp`: pristine h-AlN unit-cell/small-cell data file used by the representative script.
- `pristine/AlN_width_13.7nm.lmp`: one canonical copy of the structure used in all supplied thermostat-bias sensitivity runs (15 source copies were byte-identical).
- `doped/`: supplied equilibrated B/Ga substitution structures.
- `pores/`: supplied representative 3%, 8%, and 12% pore structures.
- `vdos/`: supplied pristine and nominal-10% B/Ga structures for VDOS simulations.

See `docs/REPRODUCIBILITY_NOTES.md` for realized-composition and pore/VDOS correspondence cautions.
