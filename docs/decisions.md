# Owner decisions

Decisions after the plan (docs/plan.md, section 9). Newest first.

## 24 Nov 2026

17. **Publication approval.** Publish the v1.0 limited preview as version 1.0.0b1 at github.com/erenextlabs/spinivo, built under decision 15 after two privacy scans of the public candidate. The pre-release number matches the open validation gates. Package index publication is a separate decision.

<!-- -->

16. **Fixture metadata.** Keep committed JEOL and Bruker fixtures byte-identical to their public deposits, including header names, account names, host names, instrument fields and paths. The deposits already publish these fields. Byte identity preserves the fixture checksums and keeps the credits accurate.

<!-- -->

15. **Fresh public history.** Start the public repository with a single initial commit containing only the files cleared by the export check, their required credits and the licence. Use Erenext LLC and the decision 12 contact for both author and committer. Include no commit trailers. Preserve the private development repository and its history without a history rewrite or force push.

<!-- -->

14. **Baseline-shape term is a diagnostic only.** Report the baseline-shape uncertainty of the area ratio in JSON and text, but keep it out of the checked-purity budget; the missing baseline-shape evidence stays an open gate of the limited preview (decision 7). Candidates were fixed before any accuracy result: the production baseline, the median of accepted quiet anchors, and a straight line through the nearest accepted anchor on each side; statistic: largest relative departure from nominal divided by √3; same-tube ratios per candidate, external-calibration tubes combined in quadrature. Measured in memory: FID06 pilot 0.13 to 0.34 %, which does not cover the −1.74 % and −1.12 % workbook deviations; the anchor-median candidate missed the FID02 anchors beside H146 by 15 and 56 noise SD and inflated the 18 EC development pairs to 4 to 15 %.

<!-- -->

13. **Holdout acceptance gates (frozen before any holdout opens; replaces the common 0.5 point target of decision 6).** Accuracy is the median absolute error against an independent reference, in mass-fraction percentage points, with repeats averaged within each preparation and each preparation weighted equally.
    - Targets: internal standard ≤0.50 points; external calibration ≤0.75 points (development EC mean absolute error 0.59 points; maleic acid tube transfer alone about 0.6 points over 18 pairs).
    - Independent preparations: at least 8 scored per method. Shared weighings, stocks, calibrants or reference lots count as one; repeats, signals and delay strata do not add. A one-sided exact sign test on the frozen roster must also pass (7 of 8 below target gives p = 0.035); ties and abstentions count as failures.
    - Fraction assessed: at least 90 % of all frozen eligible preparations per method receive a checked purity; abstentions and preparations without runs stay in the denominator.
    - Uncertainty agreement: at least 90 % of scored preparations have |error| ≤ k·u(error), k = 2 only where its coverage basis is justified. Missing reference coverage, material uncertainty or required covariance leaves this gate open. Every miss is reported.
    - Exclusions: eligibility screening only before the roster is frozen (exposed or possibly shared development material, out-of-scope experiments, no independent assigned value; aliases and repeats merged). After freezing, flags, missing inputs, processing failures and large errors stay on the roster.
    - A holdout below these counts follows decision 7 (limited preview with the gate open). A holdout search on 24 Nov 2026 (Zenodo, Harvard Dataverse, Figshare, web) found no qualified deposit.

<!-- -->

12. **Package author and contact.** Author "Erenext LLC", contact spinivo@erenext.com, project URLs on github.com/erenextlabs/spinivo. The package stays private with a development version until publication is approved.

## 23 Nov 2026

1. **Entered areas.** Purity from entered areas alone is computed but marked "not assessed". It never reports a checked pass, because it carries no acquisition or regional evidence.
2. **External transfer.** A ratio of areas from different tubes needs a traceable pair record (or a documented analyst attestation) tied to both runs and preparations, plus separate receiver gain and scale evidence. If either is missing, no purity is given.
3. **Integral checker evidence.** Only predeclared, sourced evidence (impurity band, blank, component) with alignment and scale. Without it, the checker abstains ("cannot be assessed").
4. **Acquisition evidence.** A header calibration counts only when documented as applying to the run; otherwise the user supplies a calibration record. D1 and 5 × T1 are shown as information; the recovery model decides.
5. **Holdout.** Possibly shared weighings and already seen tubes stay in development data. The final test needs another qualified public deposit or a publicly released purpose-prepared set.
6. **Acceptance target.** 0.5 percentage points median purity error for each method, provisional. Minimum independent preparations, exclusions and coverage are fixed before the holdout opens.
7. **Open gates.** If evidence stays missing, release as a limited preview that states each open gate. Nothing unverified is called qualified.
8. **Interface.** Versioned JSON plus text. CLI exit 0 for a valid report (abstention included), 1 for unsupported input, 2 for usage, 3 for a bad record.
9. **Publication content.** Only rights-cleared files: drop the vendor pulse program text, keep minimal fixtures, author and contact as approved by the owner. History cleanup and publication each need a separate approval.
10. **Baseline lock.** Adopt the locked-anchor spike test only under QUANTITATIVE_SETTINGS. Measured on 59 real regions: 10 flags before, 1 after (9 anchor-flip false alarms removed, none added), nominal outputs bit-identical. Under SETTINGS it would add 5 flags caused by the spike tail, so SETTINGS stays unchanged.
11. **Single gate.** A checked purity pass comes only from one call that reads the file, measures, checks, links and calculates. Direct calculation from entered or edited records reports purity as "calculation only" and never passes. Chosen after three reviews kept finding direct-call paths with stale or duplicated evidence.
