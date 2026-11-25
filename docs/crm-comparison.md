# Laboratory CRM comparison

Compare one checked qNMR assay with the certified mass fraction for that CRM lot.

```sh
spinivo compare comparison.json
spinivo compare comparison.json --format json
```

In Python, call `spinivo.crm_comparison.compare(path)`. Pass the report to `render_text`
or `render_json` in that module.

## Input

```json
{
  "schema": "spinivo.crm-comparison/1",
  "comparison_id": "prep-01",
  "assay_record": "assays/prep-01.json",
  "crm": {"id": "CRM-123", "lot": "A04", "certificate_source": "CRM-123 lot A04 certificate, revision 2",
          "mass_fraction_pct": 99.5, "standard_u_pp": 0.15},
  "independence": {"independent_of_calibration": true, "basis": "CRM assigned independently of the calibrant"},
  "comparison_k": {"value": 2, "basis": "Coverage factor specified in laboratory protocol 7"},
  "lab_limit_pp": {"value": 0.5, "basis": "Protocol 7 absolute-error limit fixed before measurement"}
}
```

`assay_record` is relative to the comparison file's folder. The [assay record](assay-record.md)
uses its own folder for its paths. CRM mass fraction is in percent, greater than zero and at most 100.
Uncertainties and the lab limit are in mass-fraction percentage points (pp). For expanded uncertainty,
supply `expanded_u_pp` and a positive certificate `coverage_factor` in place of `standard_u_pp`.
Spinivo divides the expanded uncertainty by that factor. Use nonnegative uncertainties and
positive values for comparison `k` and the lab limit.

## Report and verdicts

The JSON report uses `spinivo.crm-comparison-report/1` and includes `spinivo_version`. Both formats
keep the assay ID, method, status and reasons. A checked assay pass supplies the measured mass
fraction and combined standard uncertainty. Signed error is measured minus certified mass
fraction. With documented independence, the two standard uncertainties combine in quadrature
to give `u(error)`. The ratio is unavailable when `k u(error)` is zero.

- **Uncertainty agreement** passes when absolute error is at most `k u(error)`. It requires
  both standard uncertainties, comparison `k` and its basis.
- **Laboratory limit** passes when absolute error is at most the lab limit. It requires the
  predeclared limit and its basis.

Both verdicts treat absolute error as equal to the threshold when `math.isclose` returns true with
`rel_tol=1e-12` and `abs_tol=1e-12`. Both require a checked assay pass, the CRM mass fraction and a declaration
of independence with its basis. Missing or null evidence withholds the affected verdict with a reason.
The comparison and CRM IDs, lot and certificate source document provenance.

Exit status 0 means a report was produced, including failed or withheld verdicts. Unsupported
spectra exit 1, command usage errors exit 2, and bad or unreadable records exit 3. The library
raises `RecordError` for malformed data or conflicting uncertainty forms.

Neither verdict qualifies a method or releases a batch. [Decision 20](decisions.md)
