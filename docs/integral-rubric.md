# Step 4 integral labelling rubric

Pilot rubric v1. Label one assigned signal in one run at fixed bounds with a declared ¹³C satellite convention. A change of bounds creates a separate case. Record labels and their evidence before testing any checker.

Materiality is provisionally **more than 1% of the intended regional target area**. Move each edge independently inward and outward by **max(SW/N, 5% of window width)**, where N is the number of acquired complex points. The step size and the 1% bounds rule on narrow windows remain provisional pilot settings. Narrow windows omit Lorentzian wings by design. The pair truncation rule may later replace the absolute bounds rule for those wings.

One supported material defect makes a region **unreliable**. Confirmed one-sided satellite capture breaks the declared convention at any area loss. Label it **unreliable** with cause `satellites`, expect checker status **flagged**, and exclude it from purity calculation. A **reliable** label requires support for every applicable cause. Otherwise label it **indeterminate**. Preserve evidence gaps alongside supported defects. Purity certificates, workbook ratios and checker votes are insufficient regional evidence.

| Cause | Reliable evidence | Unreliable evidence | Indeterminate evidence |
|---|---|---|---|
| Neighbour or tail | A component reference or matched control bounds the foreign contribution. Assigned multiplet lines belong to the target. | An independently assigned foreign line or tail contributes material area inside the bounds. | The shoulder, extra maximum or spike difference lacks a regional assignment or contribution estimate. |
| Phase | A correctly phased reference or a known injected phase bounds area error after processing. | Residual local phase causes a material area error against the reference. | Only a global fit score, imaginary intensity or unexplained asymmetry is available. |
| Baseline | A known baseline, blank or justified range bounds area bias. | A known residual causes material bias, or two defensible placements differ materially. Record the latter as sensitivity. | Bracketing anchors or a justified baseline comparison are missing. |
| Bounds | Plausible edge moves preserve the target area within tolerance. | Any admissible move changes area by more than 1%, or component truth shows material omission. | The grid ends before a trial edge, or the intended outer component is unknown. |
| ¹³C satellites | Decoupling evidence or assigned satellite areas support the assay convention. | One-sided capture breaks the convention. Missing required target satellites or included foreign satellites cause material error. | The parent, coupling, decoupling effectiveness or convention is unresolved. |
| Low S/N | A high S/N reference or independent noise replicates support repeatability. | Replicate area scatter exceeds tolerance at the measured noise level. | A qualified noise window or independent repeatability evidence is missing. |
| Window versus linewidth | Component references show adequate capture, including the outer multiplet lines. | Supported shape and linewidth establish material truncation. | Width is unresolved or one singlet width cannot describe the envelope. |

Process all declared regions of a run together. The 22 Nov 2026 historical profile uses zero fill 2, LB 1 Hz, lineshape phase and automatic baseline correction. The quantitative profile keeps zero fill 2, fits lineshape phase on an LB 1 Hz copy and integrates at LB 0 Hz. It bridges nearby baseline anchors. Use the signed real integral in intensity·Hz as the nominal area. Relative diagnostics use its absolute magnitude. Invalid or nonpositive area disables relative metrics. Keep any independently supported defects.

For the constant baseline test, exclude blocks touching a declared region. When phase is fitted directly on the integration trace, exclude block centers within 40 fitted FWHM of a peak and subtract the fitted absorptive Lorentzian wings from quiet block levels. Estimate noise as the MAD of differences between adjacent eligible blocks divided by √2. Use the median residual level if at least 11 adjacent pairs remain and the residual MAD is at most twice that noise estimate. Otherwise use the smooth fit to the original quiet block levels. Record the nominal and shifted branch decisions, block size, eligible block and pair counts, median residual, both MAD-based statistics and the fallback reason. Repeat the baseline calculation with the grid of blocks shifted by half a block for the placement comparison.

| Cause | Proposed pilot observation |
|---|---|
| Neighbour | Use assigned positions or an aligned, scaled control to estimate foreign area inside the bounds, including external tails. An unexpected resolved maximum above 5 noise SD triggers review. |
| Phase | Compare the applied ramp with a qualified local lineshape fit or phase reference. Wrap the residual to ±180°. Rotate a diagnostic copy and compare the area with the real baseline correction accounted for. |
| Baseline | Compare the nominal and half-block areas and record both anchor distances. On the phased trace before baseline correction, add a 2 Hz FWHM absorptive Lorentzian in separate trials at 25%, 50% and 75% of the window. Scale each to 10% of the nominal in-window area and refit with all declared regions and the original phase held fixed. Report each recovery `(trial area − nominal area) / added area` and the position score `P = abs(R25 − R75)`. For smooth fits, bridge the fitted levels at the nearest accepted anchors on both sides of each declared region without a distance limit. The placement score is `C = abs(A_straight / A_nominal − 1)`. Constant and linear fits give zero placement change. Flag `baseline_sensitive` when `max(P, C) > 0.017`. The pilot limit is the rounded July FID02 upper Tukey fence and awaits independent validation. Record a reason for any failed trial or missing required bracket. |
| Bounds | Integrate all four edge trials on the nominal processed trace. Record each signed change `(trial − nominal) / abs(nominal)`. Flag any admissible trial with magnitude above 1%. An off-grid or collapsed trial leaves incomplete evidence. |
| ¹³C satellites | With a sourced coupling or matched decoupled control, inspect both bands at f₀ ± J_CH/2. Record each as inside, clipped or outside and estimate the error under the declared convention. Flag one-sided capture even below 1%. |
| S/N | Reuse `check_snr` with peak height, sample noise SD and a declared quiet window spanning at least 32 acquired points. The proposed screening minimum is 150. |
| Capture | For an isolated Lorentzian with effective FWHM L, the fraction in [l, h] is [atan(2(h−f₀)/L) − atan(2(l−f₀)/L)] / π. A centered full window of 64 × FWHM captures 99.0%. Weight fractions by assigned areas for a supported multiplet. Use an independent width or a qualified fit on the integration trace. `measure_region` supplies FWHM only at LB zero. |

`QUANTITATIVE_SETTINGS` uses threshold version `pilot-3`. Each spike refit holds the nominal quiet blocks and baseline branch fixed, so bridge endpoints stay fixed. The refit still uses the trial's block levels. Record the block centers each trial would add or remove without the lock as observations only. The nominal fit and the half-block fit use their own anchors. The 22 Nov 2026 historical profile uses `pilot-2` with the spike refit lock off.

Report absolute smooth-tail truncation loss per region as `1 − capture`, with its assumptions and uncertainty. Flag smooth-tail truncation only for an analyte–standard pair with supported `abs(capture_analyte / capture_standard − 1) > 0.005`. Apply the provisional screen to both calibration methods. Unresolved pair evidence gives **cannot be assessed**. Regional overlap, clipped assigned lines and satellite flags retain their own evidence. Pair status leaves the independent fixture labels intact.

Each future observation records run, signal, preparation, bounds, assignment source, satellite convention, processing and threshold versions, nominal area, metric, unit, limit, evidence source and reason. A supported flag takes precedence over missing evidence. With no flag, a required evidence gap gives **cannot be assessed**. All applicable checks must pass for a region to be clean and proceed to purity calculation.

## Open questions and provisional jury answers

| Open question | Provisional answer |
|---|---|
| What regional error is material, and which edge moves represent repeat placement? | Keep 1% and the stated edge moves for the pilot. Two of three jurors preferred 0.5% materiality. Record the actual changes. |
| What distinguishes a foreign tail from an assigned multiplet line? | Require independent assignment. An aligned and normalized spike supports added overlap. Quantify the native interferent for the original region. |
| What ¹³C satellite convention applies in coupled and decoupled assays? | Handle both target satellites symmetrically. Symmetric exclusion needs a validated matched convention or correction. Effective decoupling gives a merged resonance. |
| What local evidence supports adequate phase, baseline and linewidth? | Use supported local rephasing, quiet flanks and assigned outer lines. Seek a matched blank, component control or reprocessed FID when evidence is insufficient. |
| When does agreement among signals of one molecule support a label? | A range divided by mean at or below 0.5% for assigned areas per proton supports consistency only. Shared bias can persist, so agreement alone cannot establish a reliable label. |

The blind three-provider LLM jury reported these provisional answers on 22 Nov 2026. Numerical disagreements remained on edge moves, narrow-window bounds, residual satellite limits and the cross-signal cutoff.

## Frozen pilot evidence set

[The manifest](../tests/fixtures/integral_cases.json) contains 13 base synthetic regions, one separately labelled synthetic variant, 5 real regions and 3 candidate pairs with zero scored regions. Synthetic bounds use carrier-relative Hz. Real bounds use ppm.

| Kind | Reliable | Unreliable | Indeterminate |
|---|---:|---:|---:|
| Synthetic | 3 | 8 | 2 |
| Synthetic variant | 0 | 1 | 0 |
| Real | 0 | 3 | 2 |

Synthetic truth is the analytic target area, including target satellites. The historical profile uses 1 Hz broadening and the quantitative profile uses intrinsic linewidths. Use ±100 Hz for deliberate omissions and the declared bounds otherwise. At nominal S/N 3000, reliable controls must reproduce truth within 0.5% at every seed and have area CV at most 0.3%. Indeterminate cases require area CV at most 1%. Both groups require absolute half-block change at most 1%. Edge causes use the nominal trace at the default seed. The separate variant bisects one satellite on the parent trace. The S/N 300 case requires positive areas at all eight seeds and CV above 1% with and without baseline correction.

Real labels are pilot observations. Their intended response follows the recorded decoupled acquisitions. Regional decoupling effectiveness remains unverified. In the 22 Nov 2026 historical profile (LB 1, `pilot-2`), both Bruker windows exceeded 1% on an admissible edge move. At 1.15–1.22 ppm, midpoint recovery was 103.96%, P was 0.0035 percentage points and C was 3.1641%, so placement triggered `baseline_sensitive`. At 7.64–7.69 ppm, midpoint recovery was 97.91%, P was 0.1630 percentage points and C was 0.0100%. Both scores met the 1.7% limit in that window, and its half-block area change exceeded 1%. The Bruker windows are processing controls outside independent accuracy scoring. The widened Nishizaki FID05 window contains a distinct non-target band. The real set has zero independently established reliable regions.

For DVN/91YV9Q, qualify the three unspiked and GEspiked pairs by aligning spectra, normalizing to an unaffected same-tube signal and verifying dilution. Promotion requires both experiment paths, independent assignment, a predeclared target window, difference area and uncertainty. An unchanged control window supports the spike comparison. Quantitative workbook deviations are +1.0848% for FID02 I149/H146, −1.7390% for FID06 H5/F4 and −1.1197% for FID06 I5/F4. Keep all three as ratio investigations without regional labels. The 22 Nov 2026 historical LB 1 deviations for those pairs were −1.0334%, −0.8288% and −1.1967%.

## 24 Nov 2026 quantitative retest

Run `python tools/compare_references.py --checks --quantitative DATA_ROOT` for the real windows. The tool skips the widened FID05 window when its extracted data is absent. Omit `--quantitative` to reproduce the historical checks. Run `python tools/compare_references.py DATA_ROOT` for the quantitative workbook ratios above.

| Real window | Frozen label | Quantitative status and flag causes | Midpoint | P (pp) | C | Half-block change | Largest edge move |
|---|---|---|---:|---:|---:|---:|---:|
| FID02 DMSO2 | Indeterminate | Cannot be assessed | 98.0271% | 0.0464 | 0.0000% | −0.0078% | −0.0920% |
| FID02 MA | Indeterminate | Cannot be assessed | 98.2895% | 0.2905 | 0.0000% | +0.3038% | −0.0952% |
| Bruker 1.15–1.22 ppm | Unreliable | Flagged for baseline, bounds and phase | 99.9542% | 0.0026 | 3.9662% | +1.4595% | −2.8742% |
| Bruker 7.64–7.69 ppm | Unreliable | Flagged for phase | 98.2076% | 0.8089 | 0.0000% | −0.2364% | −0.7973% |
| FID05 widened DMSO2 | Unreliable | Cannot be assessed | 99.4750% | 0.1151 | 0.0000% | +0.1351% | +0.0767% |

The [fixture tests](../tests/test_integral_fixtures.py) reproduce synthetic precision of 7/7, recall of 7/8 and assessment of 11/13. The historical profile scored 6/7, 6/8 and 11/13. All reliable controls meet the 0.5% area tolerance at every seed against both analytic references. Their area CV is at most 0.3% and absolute half-block change is at most 1%. `local_baseline` passes despite a 3.007% hump. The separately scored satellite variant still flags satellites.

For owner review, the Bruker 7.64–7.69 ppm window lost its frozen baseline and bounds causes. Phase is new relative to its frozen causes and is now its sole flag. `narrow_singlet` lost its bounds cause with a largest edge move of −0.981%. All labels and their recorded causes remain frozen.

The regional checker’s sole flag for `narrow_singlet` is capture (15.596% loss, 1% limit), while the frozen rubric flags smooth-tail truncation only for an analyte–standard pair, a difference that remains unresolved.
