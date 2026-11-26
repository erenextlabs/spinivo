# Spinivo CLI user guide

## 1. What you need

Use Python 3.11 or newer and a checkout of this repository. Run the examples from the repository root.
The checkout includes the raw data and assay records used here. [Fixture sources and licences](../tests/fixtures/README.md)

For your own assay, start with a conventional 1D proton FID and the laboratory records for that preparation.
Spinivo reads JEOL JDF files and Bruker experiment folders with `fid` and `acqus`.

## 2. Install

Install Spinivo into your Python environment from the repository root.

```sh
python -m pip install .
```

You can run `python -m spinivo` in place of `spinivo` in the examples below.

## 3. Inspect the JEOL fixture

Use `inspect` to read acquisition metadata from the committed Nishizaki FID.

```sh
spinivo inspect tests/fixtures/jeol/nishizaki_fid02_ma_dmso2_d20s_ns32.jdf
```

This is the complete metadata output (exit status `0`).

```text
vendor: jeol
source: tests/fixtures/jeol/nishizaki_fid02_ma_dmso2_d20s_ns32.jdf
nucleus: 1H (x_domain)
observe_frequency_mhz: 399.78219837825 (x_freq)
spectral_width_hz: 7494.00479616307 (x_sweep)
carrier_ppm: 5.0 (x_offset)
complex_points: 65536 (x_points)
acquisition_time_s: 8.74512384 (x_points, x_sweep)
relaxation_delay_s: 20.0 (relaxation_delay)
pulse_width_us: 7.272 (x_pulse)
flip_angle_deg: 90.0 (x_angle)
pulse_90_us: 7.272 (x_90_width)
scans: 32 (scans)
dummy_scans: 4 (x_prescans)
receiver_gain: 30.0 (recvr_gain)
temperature_k: 298.15 (temp_get)
pulse_program: single_pulse_dec (experiment)
intensity_note: Stored JEOL values. No additional intensity scaling or scan normalization applied.
missing: none
intensity_scale: 1.0
storage_scale: JEOL
```

Find the nucleus, relaxation delay, pulse width and scan count. Parentheses give the header source where available.
Use `assess` to check those settings against assay evidence such as T1 and pulse calibration.

## 4. Open the assay record

The record describes a same-tube assay. Maleic acid (`MA`) is the analyte and dimethyl sulfone (`DMSO2`) is the internal standard.
Open the [JEOL assay record](../tests/fixtures/records/jeol_internal_standard.json).

These are the record's run fields, copied verbatim.

```json
  "run": {"run_id": "FID02-Pattern1-0619", "file": "../jeol/nishizaki_fid02_ma_dmso2_d20s_ns32.jdf", "preparation_id": "FID02-Pattern1-tube", "noise_bounds_ppm": [1.55, 1.85],
          "reference": {"shift_ppm": 0, "source": "JDF x_offset 5 ppm"}, "recovery": {"calibrated_90_us": 7.272, "calibration_source": "Illustrative: JDF x_90_width", "extra_recovery_s": 0}},
```

The run names the FID, the preparation, the noise window, the ppm reference and the recovery inputs.
Relative `file` paths start at the record folder. Enter signal and noise bounds on the axis after the declared reference shift.

These are the record's signal fields, copied verbatim.

```json
  "signals": {"analyte": [{"name": "MA", "bounds_ppm": [6.17, 6.28], "evidence": null, "t1_s": 2.0, "t1_source": "Illustrative T1"}],
              "standard": [{"name": "DMSO2", "bounds_ppm": [2.90, 3.01], "evidence": null, "t1_s": 2.0, "t1_source": "Illustrative T1"}]},
```

Both signal regions have bounds and T1 values. Both have `evidence: null` for component and satellite assignments.
The record also gives proton counts, molar masses, weighed masses, standard purity and uncertainty inputs.
The weighed masses follow the original FID file name.
The T1 values, pulse calibration source, molar mass sources and standard purity source are illustrative.
All uncertainties, the weighing correlation and the coverage factor `k` are also illustrative.
For a laboratory assay, use applicable lab records for each of these inputs. [Assay record fields and units](assay-record.md)

## 5. Assess and read the report

Run `assess` with the record path.

```sh
spinivo assess tests/fixtures/records/jeol_internal_standard.json
```

This excerpt starts at the report heading and ends after `Checks` (exit status `0`).

```text
FID02-MA-DMSO2 (internal_standard, internal_standard mode): cannot be assessed
  spinivo.assessment/1, spinivo 1.2.0
  Checked purity: cannot be assessed
  Baseline shape (baseline-shape/1, diagnostic): computed, 6.5074131 % relative standard uncertainty
  Blocking reasons
    MA analyte region 1: neighbour: Unexpected peak at 6.22521 ppm. capture: Missing target lines. capture: Missing target source. satellites: Missing satellite convention. satellites: Missing satellite source
    DMSO2 standard region 1: neighbour: Unexpected peak at 3.00981 ppm. neighbour: Unexpected peak at 2.95603 ppm. capture: Missing target lines. capture: Missing target source. satellites: Missing satellite convention. satellites: Missing satellite source
    MA/DMSO2 capture_pair: Missing or unusable analyte capture. Missing or unusable standard capture
  Checks
    MA analyte region 1: cannot be assessed: neighbour: Unexpected peak at 6.22521 ppm. capture: Missing target lines. capture: Missing target source. satellites: Missing satellite convention. satellites: Missing satellite source
    DMSO2 standard region 1: cannot be assessed: neighbour: Unexpected peak at 3.00981 ppm. neighbour: Unexpected peak at 2.95603 ppm. capture: Missing target lines. capture: Missing target source. satellites: Missing satellite convention. satellites: Missing satellite source
    MA recovery: pass: Meets the pilot limit for the steady-state bound
    MA snr: pass: Meets the pilot limit
    MA resolution: pass: Meets the pilot limit
    DMSO2 recovery: pass: Meets the pilot limit for the steady-state bound
    DMSO2 snr: pass: Meets the pilot limit
    DMSO2 resolution: pass: Meets the pilot limit
    MA/DMSO2 capture_pair: cannot be assessed: Missing or unusable analyte capture. Missing or unusable standard capture
```

1. Read `Checked purity` first. Spinivo withholds checked purity for this assay.
2. Read `Blocking reasons` for each region and the capture pair. Use the table below to trace each reason to evidence.
3. Read `Checks` to see which evaluations completed. Recovery, signal-to-noise ratio (`snr`) and resolution pass for both signals with the supplied illustrative inputs.
4. Scroll further for the thresholds, run settings, region observations and measured areas. Treat those areas as inputs to the checks.

The baseline-shape diagnostic estimates how the area ratio changes under alternate baseline shapes. It sits outside the checked-purity budget and all assay checks.
The separate `baseline` region check is part of the integral checks. [Baseline shape diagnostic](assay-record.md#baseline-shape-diagnostic)

Run `spinivo assess tests/fixtures/records/jeol_internal_standard.json --format html > report.html` to save the report. Open `report.html` in a browser to read or print it. The HTML report shows each raw file by its name and hides the record, current and home folders. Error messages can still name other folders, so read a report before you share it.

## 6. Trace the blocking reasons to evidence

For a laboratory assay, use assignments, calibration records and other evidence from your laboratory's own documented sources.
Leave unsupported fields missing. If evidence remains unavailable, keep the `cannot be assessed` result.

| Blocking reason in this report | What resolves the evidence gap |
| --- | --- |
| `neighbour: Unexpected peak` at 6.22521 ppm for MA and at 3.00981 and 2.95603 ppm for DMSO2 | Review each peak on the referenced spectrum against the lab's documented assignments. Record confirmed target components in `evidence.target_lines` with `evidence.target_source`. For confirmed foreign components use `evidence.foreign_lines` with `evidence.foreign_source`. Documented satellite evidence can account for an assigned satellite. |
| `capture: Missing target lines` and `capture: Missing target source` for both regions | Supply `evidence.target_lines` for each signal with its documented component positions, areas and widths (`ppm`, `area`, `fwhm_hz`). Cite the assignment or component measurement in `evidence.target_source`. Spinivo uses these inputs and `bounds_ppm` to estimate target area outside the window. |
| `satellites: Missing satellite convention` and `satellites: Missing satellite source` for both regions | Record the lab's supported treatment in `evidence.satellite_convention` and cite its basis in `evidence.satellite_source`. For `include` or `exclude`, also supply the documented `evidence.coupling_hz` and `evidence.satellite_fraction`. Use `absent` only when the evidence supports absence. |
| `capture_pair: Missing or unusable analyte capture` and `Missing or unusable standard capture` | Resolve the target capture evidence for both regions. Spinivo calculates the pair check from both capture estimates and checks their mismatch against the reported limit. |

Run `assess` again when the laboratory can supply the evidence. Review any flags that remain after Spinivo evaluates it.

## 7. Use calculation-only mode for entered areas

The FID06 calculation record uses workbook areas and an illustrative uncertainty budget. Run `assess` on that record to reproduce the arithmetic.

```sh
spinivo assess tests/fixtures/records/fid06_calculation_only.json
```

This excerpt shows the report from its heading through the expanded uncertainty, before the budget table (exit status `0`).

```text
FID06-INH-IC (internal_standard, calculation_only mode): cannot be assessed
  spinivo.assessment/1, spinivo 1.2.0
  Checked purity: cannot be assessed
  Blocking reasons
    calculation only: no spectral evidence
  Not assessed in calculation_only mode
    raw data, processing and referencing
    acquisition checks: recovery, snr, resolution
    integral region checks
    capture pair checks
  FID06-INH-IC (internal standard, same tube): cannot be assessed
    I_a = 3.6047681
    I_s = 1.2346185
    R = 2.9197424
    calculation only: no spectral evidence
    Purity = 99.92 %, u_c = 0.47 % (calculation only, not checked)
    U = 0.95 % with k = 2
```

`I_a` is the sum of the entered analyte areas `H5` and `I5`. `I_s` is the entered standard area `G5`. `R` is their ratio.
`u_c` is the combined standard uncertainty. `U` is the expanded uncertainty with coverage factor `k`.
Read the calculated value together with its `calculation only, not checked` label.

## 8. Assess an external calibration record

External calibration uses a sample and a separate calibrant. The example record pairs the committed JEOL and Bruker FIDs.
The pairing is synthetic and cannot pass. All values in the example record are illustrative.
Use it to see how Spinivo reports missing conditions and unsupported transfer.

```sh
spinivo assess tests/fixtures/records/external_calibration.json
```

This excerpt starts at the report heading and ends at the instrument transfer reason, before the region blockers (exit status `0`).

```text
synthetic-external (external_calibration, external_calibration mode): cannot be assessed
  spinivo.assessment/1, spinivo 1.2.0
  Checked purity: cannot be assessed
  Baseline shape (baseline-shape/1, diagnostic): computed, 553.24564 % relative standard uncertainty
  Blocking reasons
    sample instrument is missing
    calibrant instrument is missing
    sample probe is missing
    calibrant probe is missing
    pulse_sequence differs without K evidence
    sample decoupling is missing
    calibrant decoupling is missing
    sample excitation is missing
    calibrant excitation is missing
    temperature differs without K evidence
    sample tube is missing
    calibrant tube is missing
    sample fill is missing
    calibrant fill is missing
    observe_frequency_mhz differs; transfer between instruments is not supported
```

For a laboratory assay, link both runs and preparations with a readable pair record or attestation.
Supply documented calibrant concentration or preparation inputs, receiver scaling and response transfer evidence.
`K` is the response correction after area normalisation. Use `K_covers` for differences your laboratory's evidence supports.
Use compatible runs on the same instrument. [External calibration record fields](assay-record.md)

## 9. Read a synthetic passing report

This passing example uses a generated FID from the test suite. It represents synthetic test data only.
The excerpt starts at the report heading and ends at the last `capture_pair` line under `Checks`.

```text
synthetic-internal (internal_standard, internal_standard mode): pass
  spinivo.assessment/1, spinivo 1.2.0
  Checked purity: 0.99919407 fraction
  Baseline shape (baseline-shape/1, diagnostic): computed, 0.0010091736 % relative standard uncertainty
  Checks
    H5 analyte region 1: pass: All region checks pass
    I5 analyte region 2: pass: All region checks pass
    G5 standard region 1: pass: All region checks pass
    H5 recovery: pass: Meets the pilot limit for the steady-state bound
    H5 snr: pass: Meets the pilot limit
    H5 resolution: pass: Meets the pilot limit
    I5 recovery: pass: Meets the pilot limit for the steady-state bound
    I5 snr: pass: Meets the pilot limit
    I5 resolution: pass: Meets the pilot limit
    G5 recovery: pass: Meets the pilot limit for the steady-state bound
    G5 snr: pass: Meets the pilot limit
    G5 resolution: pass: Meets the pilot limit
    H5/G5 capture_pair: pass: Meets the pilot check
    I5/G5 capture_pair: pass: Meets the pilot check
```

Here Spinivo derives the areas from the synthetic FID and reports checked purity after the supplied evidence passes the checks.

## 10. Use JSON output and exit codes

Add `--format json` to `inspect`, `assess` or `compare` to read the output in a script.
Try `spinivo assess tests/fixtures/records/jeol_internal_standard.json --format json`.
Read `status` and `reasons` in the `spinivo.assessment/1` report. When Spinivo withholds checked purity, `purity.value` is `null`.
Use the arrayed JEOL fixture to see how `inspect` reports unsupported input in JSON.

```sh
spinivo inspect tests/fixtures/jeol/nishizaki_fid01_pulse_width_array.jdf --format json
```

This is the complete JSON output on stdout (exit status `1`).

```text
{"detail": "JDF layout is small_two_d", "reason": "arrayed_or_2d", "schema": "spinivo.inspect/1", "source": "tests/fixtures/jeol/nishizaki_fid01_pulse_width_array.jdf"}
```

`inspect` also prints a diagnostic on stderr for this input.
For `assess` or `compare`, check the report status or verdict even when the command exits with `0`.

| Exit status | Meaning | Example |
| --- | --- | --- |
| `0` | The command produced a report, including when checked purity is withheld. | `spinivo assess tests/fixtures/records/fid06_calculation_only.json` |
| `1` | Unsupported input, an inspection read error or a closed output pipe. | `spinivo inspect tests/fixtures/jeol/nishizaki_fid01_pulse_width_array.jdf --format json` |
| `2` | Command usage error. | `spinivo assess` |
| `3` | Bad or unreadable assay or comparison record. | `spinivo assess tests/fixtures/README.md` |

## 11. Compare with your laboratory's CRM

Use `spinivo compare comparison.json` with your own comparison record and a checked assay of the certified reference material (CRM).
Record the lot, certificate source, certified mass fraction and uncertainty.
Give the basis for independence from calibration, the comparison coverage factor and the laboratory's predeclared error limit.
The report gives signed error in mass-fraction percentage points and separate verdicts for uncertainty agreement and the laboratory limit.
[CRM comparison format and verdicts](crm-comparison.md)

## 12. Limits

A checked pass means the documented checks passed for the supplied assay. Spinivo makes no general purity accuracy claim.
Integral checks use pilot thresholds. Budget completeness and uncertainty coverage on real samples remain unverified.
The laboratory is responsible for method validation and batch release. [Spinivo limitations](../README.md#limitations)
