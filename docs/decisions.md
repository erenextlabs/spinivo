# Owner decisions

Decisions after the plan (docs/plan.md, section 9). Newest first.

## 25 Nov 2026

22. **User guide update.** Publish the step-by-step user guide (`docs/user-guide.md`), its README link and its output test as one follow-up commit on the public `main` under decision 18, between versions. The package version stays 1.1.0 because the code is unchanged. The same checks as decision 21 apply. Owner poll, 25 Nov 2026.

21. **Version 1.1.0 release.** Publish the v1.1 work as version 1.1.0, a regular release without the preview label, as one follow-up commit on the public `main` under decision 18. Publish only after the candidate passes the test suite, the export and licence checks and two privacy scans. Decision 20 limits stay in force, including no general accuracy claim. Package index publication remains a separate decision. Owner poll, 25 Nov 2026.

20. **Self-provable release scope.** No new third-party data, measurements or external sign-off are required to release. Spinivo claims that it reads supported files, applies its documented checks, reproduces specified purity results and GUM propagation for supported inputs, and withholds checked purity when required evidence is missing. A checked pass means the documented checks passed; it does not establish purity accuracy, budget completeness or coverage on real samples. Spinivo makes no general accuracy claim for purity, and README and reports must not imply independently validated accuracy. Calculation-only values from entered areas keep decision 1.
    - Laboratory CRM comparison: a user measures their own certified reference material; Spinivo reports the signed error in percentage points, uncertainty agreement (|error| against k·u(error)) and, separately, the verdict against the laboratory's predeclared limit. Each verdict is withheld with a reason when its inputs are unsupported. Neither verdict qualifies a method or releases a batch; method accuracy belongs to each laboratory's own validation.
    - Release gates: reader and acquisition checks, calculation and budget checks against hand-worked examples, missing-evidence checks, CLI and library agreement, licence, export and privacy, and synthetic integral-flag scores on the frozen synthetic set (`tests/fixtures/integral_cases.json`, sha256 `a8d70433bd900b9a8e3280ea6edb4c64e961d203534a95b7d940ba611e63d951`, 13 base regions, `pilot-3` thresholds): precision at least 0.90, recall at least 0.80, and at least 0.90 assessed among reliable and unreliable labels; indeterminate labels are reported separately. Changing the scored set needs a new owner decision and reports old and new scores.
    - Real-region labels follow decision 3; the internal LLM jury gives development annotations, unresolved labels stay indeterminate, and no human expert is required. Real-region scores are reported, not gated.
    - Nishizaki results are published as development evidence, labelled not independent, with every miss: workbook discrepancies as area ratios, certificate comparisons in mass-fraction percentage points. Decision 2 external-transfer checks stay per assay.
    - The decision 14 baseline-shape term stays a diagnostic outside the budget; reports and README disclose the exclusion, and it does not block a release.
    - Decision 13 gates stay defined as an optional study under decision 19 and no longer block a release; they replace the release requirements of decision 5, decision 7 for accuracy gates and the last clause of decision 13. Before any decision 19 result, the scorer must separate shared-lot covariance from independent-unit counts. The published 1.0.0b1 keeps its preview status; the version label of the first release under this scope is a separate decision.

19. **Independent preparation (amends decision 13 before any holdout opens).** A preparation is independent when it has its own weighings, its own solutions and its own tube. Shared weighings, stock solutions or calibrant solutions still count as one. A shared certified material lot (analyte reference or calibrant) does not merge preparations. Its certificate uncertainty is a shared term: include it in each preparation's error uncertainty and count it once, not as independent, in any pooled result.

## 24 Nov 2026

18. **Follow-up public releases.** Publish each later version (1.1, 1.2 and so on) as one new commit on the public `main`, with the previous public commit as its parent. Build it from the export-cleared tree under the author, committer and trailer rules of decision 15, after the same privacy scans. Never push private history and never force push; a fix release between versions needs its own approval.

<!-- -->

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
