# Spinivo plan

## 1. Goal and non-goals

Spinivo reads supported 1D proton qNMR files and applies documented acquisition and integral checks. For supported inputs it reproduces specified purity calculations and GUM propagation. When evidence is missing, it says **cannot be assessed** and identifies what is needed. Out of scope are instrument control, compound identification, regulatory certification and batch release.

## 2. Users and main use cases

| User | Main use case |
|---|---|
| QC analyst | Check a run before reporting purity. Identify a setting, signal or calibration pair that needs review. |
| Academic researcher | Assess an archived experiment and see which measurements or preparation records would make it reusable. |
| Reference material maker | Reproduce a mass-fraction result and inspect its uncertainty contributions. |
| Laboratory | Compare checked purity with its own certified reference material (CRM) with the `compare` command. |

All users get the same checks through a Python command line tool or an importable library.

## 3. Scope

**Inputs.** Read JEOL `.jdf` files and Bruker experiment folders containing `fid` and `acqus`. Support conventional 1D ¹H pulse-acquire experiments, including the calibrated 90° acquisition with ¹³C decoupling used in the Nishizaki study. Identify unsupported variants explicitly.

**Calibration.** Support a weighed internal standard in the same tube and an external calibrant in a separate tube under matched acquisition conditions. Record the method explicitly. The user supplies signal assignments, integration regions and preparation facts that the instrument cannot provide.

**Excluded.** Other nuclei, 2D data, automatic peak assignment, deconvolution, general serial-dilution workflows, calibration transfer between instruments and a graphical interface.

## 4. Architecture, measurement models and output

| Module | Responsibility |
|---|---|
| Readers | Extract the complex time-domain signal, acquisition parameters, units and source fields. Preserve missing values and identify vendor intensity scaling. |
| Processing | Apply documented digital-filter, phase and baseline handling. Produce a referenced ppm axis and integrals. Record processing settings and region bounds. Preserve relative intensity across externally paired spectra. |
| Acquisition checks | Check D1 against applicable signal T1 values. Use the pulse sequence, acquisition time and calibrated flip angle to assess recovery. Measure S/N using a stated noise region. Compare acquired Hz per point with signal linewidth before zero filling. Apply the checks to both external runs. |
| Integral checks | Assess each selected region for neighbouring signals, phase or baseline distortion, and sensitivity to modest changes in bounds. Record the observation and affected region. An expected multiplet may form one valid integral. |
| Calibration and uncertainty | Apply the selected measurement model, check its prerequisites, and propagate standard uncertainties and covariance. |
| Reporting | The CLI renders the library result as text or versioned JSON with method, IDs, versions, quantities with units, evidence, statuses and reasons. |

### Internal standard

Let `a` denote analyte and `s` standard. `I` is integral area, `N` contributing proton count, `M` molar mass and `w` mass fraction. `m_a` is the mass of sample material weighed and `m_s` the mass of standard material weighed.

`w_a = (I_a/I_s) × (N_s/N_a) × (M_a/M_s) × (m_s/m_a) × w_s`

The assay record links both signals to one tube. It supplies their assignments, proton counts, molar masses, weighed masses, standard purity, T1 evidence and uncertainties. Both integrals must be usable. The shared acquisition scale cancels in the area ratio.

### External calibration

Link both runs and preparations by a traceable pair record or documented analyst attestation. Require separate receiver gain and scale evidence for each run (decision 2). Let `c_s` be the calibrant molar concentration, `V_a` the final volume of the sample solution and `J` the area on a common per-scan receiver scale. For a directly prepared solution that contains the weighed sample:

`w_a = (J_a/J_s) × (N_s/N_a) × c_s × (M_a V_a/m_a) × K`

Use a certified `c_s` or calculate it from the calibrant preparation: `c_s = m_s w_s / (M_s V_s)`, where `V_s` is its final solution volume. The assay record gives the concentration source, preparation quantities, units and uncertainties. Distinguish final solution volume from the volume loaded into the NMR tube.

`K` is calibrant response divided by sample response after intensity normalisation. This direction matters: lower sample sensitivity requires an upward correction. The assay record identifies the instrument, probe, pulse sequence, decoupling, temperature, pulse calibration, scan accumulation, receiver gain and scaling, tube geometry, fill height and acquisition times for each run.

- Normalise areas using the documented accumulation and receiver response. Account for any scaling already applied by the reader. A nominal gain setting alone does not establish a gain ratio.
- Check that the runs have compatible excitation, relaxation, probe loading, effective observed volume and instrument stability. Record the evidence for each material response correction. Use unity only where matching and validation justify it, and retain its residual uncertainty.
- A material mismatch without a validated correction blocks external purity. The report names the mismatch and preserves the checks completed on each run.

### Uncertainty and status

For each material input, report its value, standard uncertainty, source, sensitivity coefficient and contribution. Include weighing, standard purity or concentration, integral repeatability and material molar-mass uncertainty. For external calibration also include solution preparation, relative intensity scaling, response transfer and drift. Include covariance from shared balances, stock solutions or calibrant runs where applicable. Count each contribution once. Report combined standard uncertainty and, when justified, expanded uncertainty with its coverage factor.

Checks return **pass**, **flag** or **cannot be assessed**. A checked pass means the documented checks passed. Purity accuracy, budget completeness and coverage on real samples remain unverified. Baseline shape is reported in text and JSON as an area-ratio diagnostic outside the budget (decisions 14 and 20).

Withhold checked purity as **cannot be assessed** if a required integral is unreliable or required evidence is missing. Keep observed flags, measured areas and reasons. Entered areas give a calculation-only value under decisions 1 and 11.

## 5. Test data and ground truth

| Public source | Evidence and intended use |
|---|---|
| Nishizaki et al., Harvard Dataverse, `10.7910/DVN/QFY6L9`, CC0 | **Development evidence for both methods.** Anal. Chem., `10.1021/acs.analchem.0c02967`. Caffeine CRM (99.9%) against external dimethyl sulfone. JEOL 400 MHz, calibrated 90° pulse, ¹³C decoupling, D1 20 or 60 s. The 1.13 GB ZIP contains 1,233 JDF files and five workbooks. The paper reports 16 solutions from 3 weighing sessions. Solutions 1 to 8 also contain same-tube maleic acid. Map pairs to preparations. Compare eligible external pairs against the caffeine certificate by preparation and D1. [docs/data-inventory.md](data-inventory.md) |
| Harvard Dataverse, `10.7910/DVN/MW7WUP` | Bruker, CC0. 26 bare `fid` files without `acqus`, so weak for reader tests. Four analytes; no confirmed independent analyte values. |
| Zenodo, `10.5281/zenodo.17083860` | Bruker, CC BY 4.0. Complete experiment folders with `acqus`; best Bruker reader fixture. Fermentation titers, no independent purity. |
| Harvard Dataverse, `10.7910/DVN/HXRZIY` | CC0, unpublished study; external DMSO2 calibrant, line-shape fitting; access may need a request. Low priority. |
| Project-generated public fixtures | Record the fixture licence, generation settings and expected values. Known component areas and separations test integral behaviour. Hand-worked assays test both measurement models. Score real spectra separately. |

Publish the Nishizaki development comparisons with every miss. Label them as non-independent and express workbook discrepancies as relative area-ratio differences and certificate differences in mass-fraction percentage points. [docs/open-gates.md](open-gates.md)

Keep a manifest of deposit versions, paths, checksums, licences, attributions, preparation IDs, calibration links and expected-value sources. Use only public data with verified reuse terms. Keep copyrighted figures and full texts out of the repository.

**Real integral labels.** Use decision 3's predeclared, sourced impurity band, blank or component evidence with alignment and scale. The blind three-provider LLM jury and anonymised Delphi round give development annotations. Keep unresolved or unsupported labels indeterminate. Human expert review is optional. Report real-region scores and counts by cause separately from release scoring.

## 6. Validation and acceptance

Decision 20 uses reproducible checks on existing fixtures and hand-worked examples. No new third-party data, measurements or external sign-off are required. Freeze fixture versions, expected values and tolerances before evaluation. [docs/decisions.md](decisions.md), [docs/open-gates.md](open-gates.md)

| Release gate | Acceptance check |
|---|---|
| Readers and acquisition | For JEOL and Bruker, match fixture references for nucleus, point count, spectral width, D1, intensity scale, ppm direction, peaks and areas. Test recovery, S/N and resolution boundaries. Require diagnostics for a truncated file, missing `acqus` and an unsupported experiment. |
| Calculation and budget | Match hand-worked internal-standard and external-calibration results and budgets, including covariance, to stated rounding. Verify non-unity `K`, equivalence of summed and averaged scans, and abstention for unsupported gain or probe mismatch. |
| Missing evidence | Remove and restore each required input in passing cases, including T1, pulse calibration, noise region, standard mass, external pairing, calibrant concentration or preparation inputs, and a material uncertainty. Require **cannot be assessed** with a reason, withheld dependent purity and intact independent checks. Restoration must recover the result. |
| CLI and library agreement | Match statuses, values, units and reasons for both methods in text and versioned JSON. |
| Licence | Verify redistribution rights, licence text and required credits for the candidate files. |
| Export | Pass the allowlist check and ship the cleared tree under decisions 15 and 18. |
| Privacy | Pass the privacy scans of the public candidate tree and commit message. Verify original fixture checksums under decision 16. |
| Synthetic integral flags | On `tests/fixtures/integral_cases.json` (13 base regions, `pilot-3`, SHA-256 `a8d70433bd900b9a8e3280ea6edb4c64e961d203534a95b7d940ba611e63d951`), require precision ≥0.90, recall ≥0.80 and ≥0.90 assessed among reliable and unreliable labels. Report indeterminate cases separately. Count abstention on an unreliable label as a recall miss. Changing the scored set needs a new owner decision and both old and new scores. |

The decision 13 holdout study is optional under decision 19. Accuracy, fraction assessed and uncertainty agreement are unmeasured for both methods. Decision 13 retains the ≤0.50 and ≤0.75 mass-fraction percentage-point median absolute error targets for internal and external calibration respectively, as well as the counts, exclusions and scoring rules. Separate shared-lot covariance from independent preparation counts in the scorer before reporting decision 19 results. [docs/open-gates.md](open-gates.md#optional-holdout-study)

## 7. Milestones

| Step | Done when |
|---|---|
| 1. Inventory and pair the data | The manifest identifies vendor fixtures, Nishizaki development pairs, preparation links and reference gaps. |
| 2. Read both formats | JEOL and Bruker readers pass field, intensity-scale and axis checks on real fixtures and identify unsupported inputs. |
| 3. Process and check acquisition | Reference spectra reproduce their peaks and areas. Boundary and missing-field tests give the expected statuses. |
| 4. Check integrals | The labelling rubric and evidence are recorded. Clean, contaminated and indeterminate regions exercise distinct report paths. |
| 5. Calculate internal-standard purity | A complete same-tube assay reproduces its hand-worked result and GUM budget. A blocked assay gives the specific reason. |
| 6. Calculate external-calibration purity | A matched fixture pair runs end to end. Scaling, correction direction, uncertainty and mismatch tests pass. |
| 7. Qualify and release | Pass section 6 and obtain owner publication approval. |
| Laboratory CRM comparison ([docs/crm-comparison.md](crm-comparison.md)) | Compare checked purity with the laboratory's CRM certificate per preparation. Report signed error in mass-fraction percentage points, uncertainty agreement (\|error\| ≤ k·u(error)) and a separate verdict against the predeclared absolute-error limit. Give a reason for each withheld verdict. Method validation and batch release remain the laboratory's responsibility. |

## 8. Risks and mitigations

- **Unreliable integral labels.** Bias in two areas can cancel. Check each region and retain indeterminate labels.
- **False independence.** Track shared preparations and calibrants. Disclose the sparse same-tube reference evidence.
- **External response mismatch.** Test scaling, tube matching and drift. Withhold purity where transfer cannot be justified.
- **False precision.** Show sources for budget terms and the baseline-shape exclusion.

## 9. Decisions

**Owner decisions, 22 Nov 2026**

1. Serve everyone who runs qNMR. Deliver one Python CLI and importable library for JEOL and Bruker 1D ¹H data.
2. Support internal standards in the same tube and external calibrants in separate tubes under matched acquisition. Use Nishizaki as the main external accuracy test.
3. License Spinivo under GPL-3.0-or-later. Use Python 3.11 as the minimum. The BSD-3-Clause nmrglue dependency is allowed for Bruker if its fixture tests pass. Qualify the JEOL reader separately.
4. Run a pilot first. The owner approves numerical targets and minimum evidence before the holdout opens.
5. Keep the repository private through v1.0 validation and licence review. The owner controls publication.

Later decisions supersede parts of the original plan. Decisions 5, 6, 7 and 13 record the earlier accuracy and preview rules, and decisions 19, 20, 21 and 23 define the optional study, current release scope and version labels. The private period ended with the `1.0.0b1` preview (decision 17), decision 21 sets the first regular release, `1.1.0`, under decision 20, and decision 23 sets `1.2.0`. [docs/decisions.md](decisions.md)
