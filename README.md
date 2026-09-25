# The Pfa Reachability Bound

Code and paper draft for **"The Pfa Reachability Bound: How Threshold-Offset
Growth Rate Limits Achievable False-Alarm Calibration in Parametric SAR CFAR
Detectors."** 

This repository is a trimmed, reproducibility-focused extract from a larger
multi-detector CFAR hardware research project. It contains only what is
needed to reproduce this specific paper's derivations, figures, and tables.

## What's here

```
paper/                          the paper draft (.md) and an IEEE-formatted .docx
cfar_setup.m                    run first, every session — puts everything on the MATLAB path
_common/                        shared front end, metrics, and shape-inversion utilities
_comparison/                    the sweep scripts and this paper's own verification scripts
_comparison/Results/            the summary CSVs the paper's tables and figures are drawn from
CFAR_Weibull/, CFAR Lognormal/, CFAR Generalized Gamma/,
CFAR_G0/, CFAR_Burr/, CFAR K/   the six clutter-model estimators + threshold functions
```

Each detector folder holds exactly three files: `*_Params.m` (MoLC shape
estimator), `*_TLog.m` (the log-domain threshold offset `delta(shape, Pfa)`
this paper's growth-rate derivations are about), and `*_Floating.m` /
`*_Shared.m` (the standard detector-contract wrapper the sweep scripts call).
Fixed-point/RTL/hardware-generation files, legacy scripts, and demo drivers
that exist in the parent research project are deliberately not included —
they aren't needed for anything this paper reports.


## Reproducing the paper's results

```matlab
cfar_setup();                 % adds everything above to the MATLAB path

verify_growth_rates()         % Section 2.1 -- cross-checks every derived
                               % growth-rate formula against the real
                               % *CFAR_TLog.m implementations, from the
                               % actual sweep Pfa range down to Pfa=1e-40

compute_clopper_pearson()     % Section 4 -- attaches exact Clopper-Pearson
                               % 95% confidence intervals to every reported
                               % achieved-Pfa figure and prints the paper's
                               % headline table with its intervals
```

Both of the above run immediately from the CSVs already in
`_comparison/Results/` — no dataset download needed.

**Section 3 and the Section 4 addendum** (the BC-3/cross-scene comparison
and the 250-image HRSID sample-expansion check) were computed from
`comparison_summary_ssdd_full_fixed60.csv`, `comparison_summary_hrsid_full.csv`,
and `comparison_summary_hrsid_full250.csv`, all included in
`_comparison/Results/` — reading those directly reproduces every number in
those sections without needing to re-run anything.

**To regenerate those summary CSVs from raw images instead of trusting the
included ones** (not necessary to verify the paper's numbers, only to
re-derive them from scratch), you need the two source datasets, which are
not redistributed here:

- **SSDD** (SAR Ship Detection Dataset) — place at `_datasets/SSDD/JPEGImages/`
  and `_datasets/SSDD/Annotations/` (VOC-format XML).
- **HRSID** (High-Resolution SAR Images Dataset) — place at `HRSID/images/`
  and `HRSID/annotations/train_test2017.json` (COCO-format).

With those in place:

```matlab
run_comparison('NumImages', 60, 'Sli', [11 15 17 21 31 41 51 61 71 81 91 101 111 121 131 141 151], 'Tag', '_ssdd_full_fixed60');
run_comparison_hrsid('NumImages', 60, 'Sli', [11 15 17 21 31 41 51 61 71 81 91 101 111 121 131 141 151], 'Tag', '_hrsid_full');
run_comparison_hrsid('NumImages', 250, 'Sli', [11 15 17 21 31 41 51 61 71 81 91 101 111 121 131 141 151], 'Random', true, 'Seed', 42, 'Tag', '_hrsid_full250');
```

## Status

The paper draft is content-complete (all six roadmap items resolved: the
six-distribution growth-rate taxonomy, the cross-scene collision check,
Clopper-Pearson intervals, the HRSID sample expansion, and an audited
literature search); all reference entries are verified against their
source records. 
