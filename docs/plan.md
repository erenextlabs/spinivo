# Spinivo v1.0 plan

## 1. Goal and non-goals

Spinivo will assess whether a 1D proton qNMR experiment supports a reported purity. It will check acquisition settings and integral regions, calculate analyte mass fraction by either of two calibration methods, and report a GUM uncertainty budget. When evidence is missing, it will say **cannot be assessed** and identify what is needed.

Out of purpose: instrument control, compound identification, regulatory certification and release decisions.

## 2. Users and main use cases

| User | Main use case |
|---|---|
| QC analyst | Check a run before reporting purity. Identify a setting, signal or calibration pair that needs review. |
| Academic researcher | Assess an archived experiment and see which measurements or preparation records would make it reusable. |
| Reference material maker | Reproduce a mass-fraction result and inspect each contribution to its uncertainty. |

All users get the same checks through a Python command line tool or an importable library.

## 3. v1.0 scope

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
| Reporting | One library result; the CLI renders plain text or versioned JSON with method, IDs, versions, quantities with units, evidence, check outcomes and reasons. |

### Internal standard

Let `a` denote analyte and `s` standard. `I` is integral area, `N` contributing proton count, `M` molar mass and `w` mass fraction. `m_a` is the mass of sample material weighed and `m_s` the mass of standard material weighed.

`w_a = (I_a/I_s) × (N_s/N_a) × (M_a/M_s) × (m_s/m_a) × w_s`

The assay record links both signals to one tube. It supplies their assignments, proton counts, molar masses, weighed masses, standard purity, T1 evidence and uncertainties. Both integrals must be usable. The shared acquisition scale cancels in the area ratio.

### External calibration

Link the sample run to a calibrant run by file and preparation ID. Let `c_s` be the calibrant molar concentration, `V_a` the final volume of the sample solution and `J` the area on a common per-scan receiver scale. For a directly prepared solution that contains the weighed sample:

`w_a = (J_a/J_s) × (N_s/N_a) × c_s × (M_a V_a/m_a) × K`

Use a certified `c_s` or calculate it from the calibrant preparation: `c_s = m_s w_s / (M_s V_s)`, where `V_s` is its final solution volume. The assay record gives the concentration source, preparation quantities, units and uncertainties. Distinguish final solution volume from the volume loaded into the NMR tube.

`K` is calibrant response divided by sample response after intensity normalisation. This direction matters: lower sample sensitivity requires an upward correction. The assay record identifies the instrument, probe, pulse sequence, decoupling, temperature, pulse calibration, scan accumulation, receiver gain and scaling, tube geometry, fill height and acquisition times for each run.

- Normalise areas using the documented accumulation and receiver response. Account for any scaling already applied by the reader. A nominal gain setting alone does not establish a gain ratio.
- Check that the runs have compatible excitation, relaxation, probe loading, effective observed volume and instrument stability. Record the evidence for each material response correction. Use unity only where matching and validation justify it, and retain its residual uncertainty.
- A material mismatch without a validated correction blocks external purity. The report names the mismatch and preserves the checks completed on each run.

### Uncertainty and status

For each material input, report its value, standard uncertainty, source, sensitivity coefficient and contribution. Include weighing, standard purity or concentration, integral repeatability and material molar-mass uncertainty. For external calibration also include solution preparation, relative intensity scaling, response transfer and drift. Include covariance from shared balances, stock solutions or calibrant runs where applicable. Count each contribution once. Report combined standard uncertainty and, when justified, expanded uncertainty with its coverage factor.

Each check returns **pass**, **flag** or **cannot be assessed**. An observed defect retains its flag. If a required integral is unreliable, or material acquisition, calibration or uncertainty evidence is unresolved, the final purity is **cannot be assessed**. The report retains the measured areas and the reason for withholding purity.

## 5. Test data and ground truth

| Public source | Evidence and intended use |
|---|---|
| Nishizaki et al., Harvard Dataverse, `10.7910/DVN/QFY6L9`, CC0 | **Main external-calibration accuracy test.** Paper: Anal. Chem., `10.1021/acs.analchem.0c02967`. Caffeine CRM (certified 99.9%) against external dimethyl sulfone; JEOL 400 MHz, calibrated 90° pulse, ¹³C decoupling, D1 20 or 60 s. One 1.13 GB ZIP holds raw data and calculation sheets. Paper reports 16 solutions from 3 weighing sessions; solutions 1 to 8 also hold maleic acid as a same-tube internal standard. Archive holds 1,233 JDF files and five workbooks (see `data-inventory.md`). Map pairs and preparations before scoring. Reserve the caffeine certificate value for comparison. |
| Harvard Dataverse, `10.7910/DVN/MW7WUP` | Bruker, CC0. 26 bare `fid` files without `acqus`, so weak for reader tests. Four analytes; no confirmed independent analyte values. |
| Zenodo, `10.5281/zenodo.17083860` | Bruker, CC BY 4.0. Complete experiment folders with `acqus`; best Bruker reader fixture. Fermentation titers, no independent purity. |
| Harvard Dataverse, `10.7910/DVN/HXRZIY` | CC0, unpublished study; external DMSO2 calibrant, line-shape fitting; access may need a request. Low priority. |
| Project-generated public fixtures | Record the fixture licence, generation settings and expected values. Known component areas and separations test integral behaviour. Hand-worked assays test both measurement models. Score real spectra separately. |

**Internal-standard accuracy:** Nishizaki solutions 1 to 8 are the first candidate; a Bruker same-tube set with an independent reference is still missing. A calculation sheet shows reproducibility, not accuracy.

Keep a manifest of deposit versions, paths, checksums, licences, attributions, preparation IDs, calibration links and expected-value sources. Use only public data with verified reuse terms. Keep copyrighted figures and full texts out of the repository.

**Integral ground truth.** Label each region reliable, unreliable or indeterminate. Record independent evidence, such as a clean component spectrum or a known added interferent. A sample purity certificate cannot establish which region is reliable. The review method uses a blind three-provider LLM jury, an anonymised Delphi round for splits and human expert review for remaining splits. Factual labels need checkable evidence. Otherwise keep them indeterminate.

## 6. Validation and acceptance

Pilot first. Before opening the holdout, the owner approves targets, minimum independent case counts and permitted exclusions. Freeze rules, reference values and group assignments. Keep runs from one preparation together. Keep shared external calibration batches together where possible and disclose any dependency. Report case counts, exclusions and abstentions by vendor and calibration method.

| Test | Acceptance check |
|---|---|
| Readers and acquisition | For each vendor, compare nucleus, point count, spectral width, D1, ppm direction and a known peak position against independent records. Exercise both sides of rule boundaries. A pulse width without calibration cannot establish flip angle. |
| Integral flags | Report precision, recall, counts by failure cause and fraction assessed on held-out real regions. Count an abstention on a known unreliable region as a miss for recall. **Provisional** targets are precision ≥0.90, recall ≥0.80 and ≥0.90 assessed among labelled regions. Report synthetic scores separately. |
| Internal purity | Reproduce a hand-worked result and budget, including covariance, to stated rounding. Score real accuracy on eligible same-tube preparations once available. |
| External purity | Reproduce a hand-worked matched pair and budget. Test a known non-unity `K` for correction direction, equivalent summed and averaged scans for identical results, and an unsupported gain or probe mismatch for abstention. Compare eligible Nishizaki pairs against 99.9% caffeine by preparation and D1. Compare intermediate quantities with the sheets to locate discrepancies. |
| Accuracy for each method | Report signed bias, median absolute error in mass-fraction percentage points, agreement relative to combined measurement and reference uncertainty, and the fraction of eligible preparations given a result. Median absolute error targets are ≤0.50 percentage points for internal standard and ≤0.75 for external calibration; counts, coverage and exclusions are in owner decision 13 (docs/decisions.md). Record reference uncertainty and its coverage basis before judging uncertainty agreement. |
| Missing evidence | Remove one required fact at a time from passing assays. Include T1, pulse calibration, noise region, standard mass for a weighed standard, external run pairing, calibrant concentration or its preparation inputs, and a material uncertainty term. Require **cannot be assessed**, the specific reason, no dependent purity value and intact independent checks. Restoring the fact must restore the original result. |
| End to end | Match CLI and library statuses, values, units and reasons for both methods. Require specific diagnostics for a truncated file, missing `acqus` and an unsupported experiment. |

If a holdout lacks enough independent preparations or labelled regions, keep the affected release check open.

## 7. Milestones

| Step | Done when |
|---|---|
| 1. Inventory and pair the data | The manifest identifies qualified vendor fixtures, eligible Nishizaki pairs, independent preparations and gaps in same-tube reference data. |
| 2. Read both formats | JEOL and Bruker readers pass field, intensity-scale and axis checks on real fixtures and identify unsupported inputs. |
| 3. Process and check acquisition | Reference spectra reproduce their peaks and areas. Boundary and missing-field tests give the expected statuses. |
| 4. Check integrals | The labelling rubric and evidence are recorded. Clean, contaminated and indeterminate regions exercise distinct report paths. |
| 5. Calculate internal-standard purity | A complete same-tube assay reproduces its hand-worked result and GUM budget. A blocked assay gives the specific reason. |
| 6. Calculate external-calibration purity | A verified Nishizaki pilot pair runs end to end. Scaling, correction direction, uncertainty and mismatch tests pass. |
| 7. Qualify and release v1.0 | The owner approves targets before holdout access. Both methods meet their frozen gates, CLI and library agree, and licence review is complete. The owner then authorises publication. |

## 8. Risks and mitigations

- **Unreliable integral labels.** Bias in two areas can cancel in final purity. Validate individual regions and retain indeterminate labels.
- **False independence.** Many spectra may share one preparation or calibrant. Track those links and report independent counts.
- **External response mismatch.** Test intensity handling, tube matching and drift on paired runs. Withhold purity where transfer cannot be justified.
- **Sparse same-tube evidence.** Check candidate preparations early. Keep the internal-standard accuracy gate open until eligible public cases exist.
- **False precision.** Show sources for material budget terms and test missing uncertainties explicitly.

## 9. Decisions

**Owner decisions, 22 Nov 2026**

1. Serve everyone who runs qNMR. Deliver one Python CLI and importable library for JEOL and Bruker 1D ¹H data.
2. Support internal standards in the same tube and external calibrants in separate tubes under matched acquisition. Use Nishizaki as the main external accuracy test.
3. License Spinivo under GPL-3.0-or-later. Use Python 3.11 as the minimum. The BSD-3-Clause nmrglue dependency is allowed for Bruker if its fixture tests pass. Qualify the JEOL reader separately.
4. Run a pilot first. The owner approves numerical targets and minimum evidence before the holdout opens.
5. Keep the repository private through v1.0 validation and licence review. The owner controls publication.

Later owner decisions: docs/decisions.md.

### Still open

| Decision | Options and recommendation |
|---|---|
| Final acceptance values and case counts | **Decided 24 Nov 2026:** separate gates for each method, owner decision 13 in docs/decisions.md. |
| Source of internal-standard accuracy cases | **Qualify an existing public deposit first**, or seek a purpose-prepared assay released publicly with its preparation record and independent reference. Existing candidates can be assessed as soon as deposit verification finishes. |
