# Literature validation benchmarks for pristine monolayer h-AlN

Present 300 K NEMD finite-size extrapolations in the supplied workbooks are approximately **116.97 W m^-1 K^-1 (armchair)** and **117.17 W m^-1 K^-1 (zigzag)**. Literature values are strongly method-dependent, so the table is intended as a methodological comparison rather than a claim that any one published value is uniquely correct.

| Reference | Method | Reported/representative room-temperature h-AlN kappa | Use in discussion |
|---|---|---:|---|
| Kocabaş et al. (2023), DOI 10.1039/D3NR00399J | DFT/GAP/HIPHIVE force constants + BTE | ~108-118 W m^-1 K^-1 across several routes | Closest multi-method benchmark to the present ~117 W m^-1 K^-1 result. |
| Zhang et al. (2025), DOI 10.1016/j.commatsci.2025.113906 | First-principles phonon transport / BTE-Wigner framework | approximately ~119 W m^-1 K^-1 for the relevant reported treatment in the project validation notes | Independent first-principles comparison. |
| Karaaslan et al. (2020), DOI 10.1103/PhysRevApplied.13.034027 | Classical MD; transferable Tersoff-type group-III-nitride potential | Potential/method benchmark rather than the sole target value | Primary provenance for the potential family and classical-MD thermal-transport methodology. |
| Karaaslan et al. (2021), DOI 10.1063/5.0051975 | Classical MD; defect thermal transport; GK/NEMD comparisons | Methodological/defect benchmark | Useful for vacancy/pore and NEMD methodology context. |
| Wang et al. (2021), DOI 10.1088/1361-6528/abd20c | First-principles BTE | 74.42 W m^-1 K^-1 for monolayer AlN | Demonstrates lower values reported with other first-principles treatments. |
| Banerjee et al. (2022), DOI 10.1039/D2CP01513G | Ab-initio iterative BTE / ShengBTE | 306.5 W m^-1 K^-1 at 300 K | Demonstrates a much higher first-principles prediction and a strong tensile-strain enhancement. |
| Duan et al. (2025), DOI 10.1016/j.commatsci.2025.113987 | DFT + iterative BTE | low room-temperature value reported in the source validation notes | Demonstrates the broad method-dependent spread in published h-AlN predictions. |

Do **not** describe higher/lower literature results as intrinsically erroneous. Differences in force constants/potentials, anharmonic treatment, size assumptions, thickness convention, BTE implementation, and other modeling choices can materially change reported 2D thermal conductivity.
