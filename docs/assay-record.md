# Assay record `spinivo.assay/1`

One JSON file describes one assessment. `spinivo.load_record(path)` reads it, `spinivo.parse_record(data, base_dir)` checks a decoded object and `spinivo.assess_record(record)` runs the route the mode names. Examples: [`tests/fixtures/records/`](../tests/fixtures/records/).

- **Bad record.** Invalid JSON, numbers outside the finite float range, unknown keys, wrong types or list lengths, unknown modes, blank or repeated signal names, and spectral runs without `file` raise `RecordError` with the field path. Component list items must be objects.
- **Missing facts.** An absent field is null and an absent list is empty. Missing evidence loads and gives "cannot be assessed" with a named reason. State `extra_recovery_s` even when it is 0. Supply `satellite_fraction` for `include` or `exclude`.
- **Gate.** The spectral modes return the assessment of `assess_internal_standard` or `assess_external_calibration` (decision 11) and accept no areas or check statuses. `calculation_only` takes entered areas for the internal-standard model and returns the direct calculation result, which is always "cannot be assessed" (decision 1).
- **Run ID.** Supply a `run_id` for every spectral run. For a direct call to `assess_internal_standard`, pass `run_id` as a keyword argument. If the ID is missing, the internal assessment reports "run_id is missing" and withholds purity.
- **Paths and pair.** Relative `file` and `pair.source` paths are bound to the absolute record folder (`base_dir`). Use `file: "."` when that folder is the experiment. The pair names both preparations and two distinct run labels that match the entered `run_id` values before they map to file paths. The external pipeline requires `pair.source` to be a readable, nonempty file outside both run folders.
- **Reference.** Declare `reference.shift_ppm` with its source for every spectral run. Processing adds the shift to the header ppm axis. Enter signal bounds, noise bounds and component positions on the resulting axis.

| Object | Fields |
|---|---|
| top level | `schema` (`"spinivo.assay/1"`, required); `mode` (`internal_standard`, `external_calibration` or `calculation_only`, required; `calculation_only` uses the internal-standard fields); `notes` (free text for sources and licences, not used) |
| | `assay_id`; `N_a`, `N_s` (proton counts of the analyte or sample and of the standard or calibrant signals); `R_relative_u`, `R_u_source`, `R_source_ids` (relative standard uncertainty of the area ratio R, its source and shared source IDs); `k` (coverage factor `{value, source}`) |
| | `quantities`: quantity objects by name. Internal standard: `M_a`, `M_s`, `m_a`, `m_s`, `w_s`. External: `M_a`, `V_a`, `m_a`, `K`, `D`, and `c_s` (certified route) or `m_s`, `w_s`, `M_s`, `V_s` (prepared route) |
| | `correlations`, `covariances`: lists of declarations between budget inputs (`R` and the quantity names) |
| | `signals`: signal lists by role, `analyte` and `standard` or `sample` and `calibrant`. Runs: `run` (internal standard, the one tube) or `sample` and `calibrant` (external) |
| | External only: `c_s_route` (`certified` or `prepared`); `K_covers` (condition names whose differences K justifies); `pair`; `scale` |
| quantity | `value`, `standard_u`, `unit`, `source`, `source_ids` (shared sources, such as one balance) |
| correlation | `inputs` (two input names), `r`, `source`, `independence_reason` (for a zero value) |
| covariance | `inputs`, `value`, `units` (both input units in input order), `source`, `independence_reason` |
| signal | `name` (nonblank, unique in the record). Spectral modes: `bounds_ppm` (two bounds), `evidence`, `t1_s`, `t1_source`. `calculation_only`: `area`, `area_source` |
| evidence | `target_lines`, `target_source`, `foreign_lines`, `foreign_source` (sourced component lines `{ppm, area, fwhm_hz}`); `satellite_convention` (`include`, `exclude` or `absent`), `coupling_hz`, `satellite_fraction`, `satellite_source` |
| run | `run_id` (run label), `preparation_id` (tube or solution), kept apart from the file. Spectral modes: `file` (JDF file or Bruker folder), `noise_bounds_ppm` (two bounds), `reference` `{shift_ppm, source}`, `recovery` `{calibrated_flip_angle_deg or calibrated_90_us, calibration_source, extra_recovery_s, extra_recovery_source}` |
| | External: `conditions` with text `area_mode`, `probe`, `pulse_sequence`, `decoupling`, `solvent`, `tube`, `instrument`, `excitation`, `fill`, `reader_scale_source` and numbers `scans`, `receiver_gain`, `temperature` (°C), `reader_scale`, `relaxation_delay_s`, `acquisition_time_s`, `spectral_width_hz`, `complex_points`, `observe_frequency_mhz` |
| pair | `kind` (`record` or `attestation`), `source` (path of the pair record), `sample_run_id`, `sample_preparation_id`, `calibrant_run_id`, `calibrant_preparation_id` |
| scale | `direction` (`"sample/calibrant"`, required when `scale` is given); `S` (a quantity); `gains`, `storage_scales` (the receiver gains and storage scales S covers, in sample, calibrant order) |

## Baseline shape diagnostic

Spectral assessment reports include `baseline_shape` and one text line. The JSON diagnostic gives
`method_version`, `status`, `reason`, `relative_standard_uncertainty_percent` and `regions`. Each region
identifies its run, role, signal and ppm bounds and gives candidate areas in intensity·Hz.
The candidates are the production baseline, the median level of all accepted quiet anchors, and
the line through the nearest accepted anchor on each side. Phase, bounds and anchors stay fixed.

For each candidate, sum region areas by role. For a same-tube assay, form the ratios from those
sums and divide the largest relative departure from the production ratio by √3. For external
calibration, apply the same statistic to each tube's summed areas and combine the two relative
terms in quadrature, assuming independent residual baseline errors across tubes.

The diagnostic is `computed` under `QUANTITATIVE_SETTINGS` when the candidate areas are usable.
Other processing settings give `not computed`. A missing anchor on either side or a nonpositive
or nonfinite area gives `cannot be assessed`. The report gives a reason and retains available
candidate areas. Unavailable and nonfinite values appear as JSON null. Calculation-only reports
omit the diagnostic. Decision 14 keeps this diagnostic outside the purity result, uncertainty
budget and all assay checks.
