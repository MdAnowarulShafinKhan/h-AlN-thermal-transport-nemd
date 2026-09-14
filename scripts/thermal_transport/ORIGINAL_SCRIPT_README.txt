FINAL UNIFIED NEMD THERMAL-TRANSPORT ANALYSIS SCRIPT
====================================================

Script:
  thermal_transport_analysis_ALL_CASES_FINAL.py

Normal command:
  python thermal_transport_analysis_ALL_CASES_FINAL.py your_file.xlsx

Default behavior:
  - 10,000 bootstrap iterations
  - 1200 dpi standalone PNG figures
  - first worksheet is analyzed
  - output: your_file_analyzed.xlsx
  - figures: your_file_analysis_figures/

Install dependencies:
  python -m pip install -r thermal_transport_ALL_CASES_requirements.txt

SUPPORTED INPUTS
----------------

A. Non-doped, one or multiple lengths

  Length (nm) | Seed 1 | Seed 2 | Seed 3 | Seed 4 | Seed 5

  One row:
    Calculates finite-length mean kappa and two-sided 95% Student-t CI.
    Does not calculate kappa_inf or lambda_eff.

  Three or more length rows:
    Calculates point statistics, reciprocal fit, kappa_inf, lambda_eff,
    10,000-iteration bootstrap CIs, residuals, shortest-length exclusion,
    and leave-one-length-out sensitivity.

B. Doped, five realizations x one seed

  Length (nm) | r1 | r2 | r3 | r4 | r5

or

  Length (nm) | R1_S1 | R2_S1 | R3_S1 | R4_S1 | R5_S1

  One row gives finite-length kappa and 95% CI across realizations.
  Multiple genuine lengths give kappa_inf and lambda_eff.

C. Doped, five realizations x three seeds

  Length (nm) | R1_S1 | R1_S2 | R1_S3 | ... | R5_S3

  The script:
    1. averages the three seeds inside each realization;
    2. uses the five realization means as the independent statistical units;
    3. calculates the Student-t 95% CI across the five realization means;
    4. uses hierarchical bootstrap resampling of realizations and seeds;
    5. performs the seed-noise audit.

  The default preselected seed for the reduced one-seed audit is S1.
  To audit S2 instead:
    python thermal_transport_analysis_ALL_CASES_FINAL.py file.xlsx --preselected-seed 2

D. Hierarchical long format

  Length (nm) | Realization | Seed | Kappa

E. Single-length or screening series

  Accepted first-column types include:
    Width, Delta, Strain, Concentration, Porosity, Temperature, Direction.

  These receive point-wise mean kappa and two-sided 95% Student-t CIs.
  Width, Delta, and two-direction data also receive automatic comparison checks.
  No kappa_inf or lambda_eff is calculated because these are not length scaling.

IMPORTANT RULES
---------------

1. The script expects individual thermal conductivities as inputs. It does not
   calculate kappa from raw heat flux and temperature-gradient files.

2. kappa_inf and lambda_eff are calculated only when the first column is a
   genuine transport length and at least three unique lengths are supplied.
   Five or more lengths are recommended and flagged more favorably.

3. Two lengths are rejected because they cannot support a defensible fit
   diagnostic.

4. For formal doped length scaling after a successful seed-noise audit, use one
   preselected seed per realization at every length. If the audit fails and all
   lengths are expanded to three seeds per realization, use the hierarchical
   format at every length.

5. Use --paired-realizations only if the realization IDs are deliberately
   matched across all lengths. Otherwise omit it.

6. Use --paired-units only if simple replicate columns represent deliberately
   matched units across all lengths. Ordinary seed replicates are normally
   unpaired, so omit it.

7. Sheet conductance:
    python thermal_transport_analysis_ALL_CASES_FINAL.py file.xlsx --thickness-nm VALUE

8. Exact realized dopant concentration:
    python thermal_transport_analysis_ALL_CASES_FINAL.py file.xlsx \
      --dopant-count N_DOPANTS --host-sites N_HOST_SITES

9. Select a worksheet:
    python thermal_transport_analysis_ALL_CASES_FINAL.py file.xlsx --sheet SheetName

10. Override a nonstandard x column or interpretation:
    --x-column "exact header"
    --analysis-kind length
    --analysis-kind screening
    --analysis-kind width
    --analysis-kind delta
    --analysis-kind direction

11. Default figures are 1200 dpi. Override when needed:
    --dpi 600

OUTPUT WORKBOOK
---------------

The workbook includes, as applicable:
  README
  Raw_Data
  Realization_Means
  Point_Statistics
  Point_Bootstrap_Summary
  Fit_Summary
  Fit_Bootstrap_Summary
  Recommended_Reporting
  Sensitivity
  Seed_Noise_Audit
  Special_Checks
  Diagnostics
  OriginPro_Data
  Bootstrap_Samples
  Figures

SCIENTIFIC LIMITS
-----------------

The script does not:
  - calculate VDOS/PDOS;
  - inspect raw temperature-profile linearity;
  - inspect cumulative thermostat-energy linearity;
  - create one project-wide comparison chart from many separate case files.

lambda_eff is an effective characteristic length from
  1/kappa(L) = a + b/L,
  kappa_inf = 1/a,
  lambda_eff = b/a.
It is not a mode-resolved phonon-MFP spectrum.
