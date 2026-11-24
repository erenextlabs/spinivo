# Data inventory

Checked 22 Nov 2026 on checksum-verified downloads (MD5 matches the deposits). Paths are relative to the local `.cache/veri/` copy, which is not in the repo.

## Nishizaki · 10.7910/DVN/QFY6L9 · CC0

1,233 JDF files and five calculation workbooks, FID02 through FID06. CAF = caffeine. MA = maleic acid. EC = external calibrant.

| Folder | JDFs | Run map and workbook links |
|---|---:|---|
| FID01_Table2 | 126 | MA + CAF, filename recipe matches solution 04. Days 1–3 × delays 1, 5, 10, 20, 30, 60, 120 s × CAF H1/H3/H7/H8, MA and solvent signal. |
| FID02_Figure2 | 816 | 5 mm tubes 01–04 and 3 mm tubes 05–08. Each tube size has Patterns 1–6 and six dates per pattern. Pair CAF 04 with DMSO2 01/02/03, and CAF 08 with DMSO2 05/06/07, within tube size, pattern and date. P1 is D20/NS32 with fixed pulse width. P2 has qH and nutation, with D1 absent from filenames. P3–6 generally have D20 qH and nutation. P5 also has 3 mm `1h_90` and `qH_auto90pw_tune` files. |
| FID03_Table3 | 72 | 5 mm, D20/NS16. Patterns 1/2/4 × N1–N3. CAF 09 is EC for CAF 10 and conductive CAF 11 within each group. P1 uses a fixed pulse, P2/4 include nutation. The sheet uses 10's area in every 11 formula while taking concentration and pulse width from 09. |
| FID04_Figure3 | 96 | 3 mm, 05–08, Patterns 7/8 × N1–N6, RG40/NS32. P7 uses `SHIMqH`, P8 `qH`. Within each group pair CAF 08 with DMSO2 05/06/07. Sheet headers and concentrations identify 05–08 despite row IDs 1–4. |
| FID05_Figure4 | 111 | Pattern 4. In each decoupling mode and N1–N3, pair DMSO2 13–16 at D60 with DMSO2 12 at D20 or D60. Solvents for 12–16 are DMSO-d6, CD3OD, acetone-d6, CDCl3 and D2O. Includes four T1 files for 13–16. Non-decoupled N3 has an apparent misplaced D2O nutation file named `12` under `Delay20sec`. |
| FID06_Table4 | 12 | D60/NS16, qH and nutation. Isoniazid + BTMSB pairs with DMSO2 13 in CD3OD. Rifampicin + DSS-d6 pairs with DMSO2 12 in DMSO-d6. The latter EC files are byte-identical to FID05 non-decoupled N1/D60. |

- **Same-tube standards and gravimetry.** FID02 confirms MA in 01–08 by concentration, integral and same-row IC formulas. FID04 reuses 05–08. FID06 adds BTMSB and DSS-d6. FID02–05 give concentrations and calculated results, with no explicit balance records or purity inputs. Filename values are consistent with encoded masses and volumes. For example, 01's 4.226 mg MA and 1.534 mg DMSO2 in 6 mL reproduce the sheet at 99.9%, as does 04's 15.130 mg CAF in 9 mL. FID06 formulas explicitly use 1.064/8.463 mg BTMSB/isoniazid and 0.861/9.817 mg DSS-d6/rifampicin, each in 1.5 mL. Its footnote gives 99.9% DMSO2 and BTMSB, 92.4% DSS-d6.
- **Three headers.** The two FID02 pilot files for 04 and 01 below are 1D proton, 65,536 points, SW 7,494.005 Hz (clipped 5,995.204 Hz), carrier 399.782 MHz, offset 5 ppm, D1 20 s, NS32, RG30, probe ID 3448 and 25 °C setpoint. Both record `single_pulse_dec` with carbon-13 decoupling. `x_pulse` is 7.36572/7.37538 µs, matching the sheet. FID05 sample 13 below gives D1 60 s, NS16 and 7.27788 µs, also matching the sheet. Acquisition time, prescans, solvent and filter fields are recoverable.
- **Eligibility.** Reader fixture yes. Acquisition checks yes for settings, delay and pulse contrasts, with nutation and T1 data to link. Internal accuracy candidate for CAF 04/08 and DMSO2 01–03/05–07 versus CRM values after preparation verification. FID06 supports IC/EC calculation reproduction. External accuracy candidate for 04/08 versus 99.9% CAF via DMSO2. FID03 tests conductivity and FID05 tests solvent, delay and decoupling effects. Integral labels open. Sheet areas support reproduction and FID05's separate-solvent spectra can help assess overlap.

### Candidate JEOL pilot pairs

Paths below are relative to `.cache/veri/nishizaki_QFY6L9/x_334_Dataverse_20200705/334_Dataverse_20200705/`. Join each directory with its two filenames.

- CAF accuracy and same-tube MA pilot, FID02 5 mm sheet rows 146/149, D20.
  `FID02_Figure2/01_5mmNMRtube/Pattern4/0724_2019/`
  - Sample `04_MAandCAF_T1A_5mm_6435u15130u9000uL_yuzo190724_RG30_D20s_qH-1-1.jdf`
  - EC `01_MAandDMSO2_T1A_5mm_4226u1534u6000uL_yuzo190724_RG30_D20s_qH-1-1.jdf`
- D60 matched-delay solvent transfer, FID05 cells J7/K7.
  `FID05_Figure4/01_13C_decoupled_Aquisition/N1/Delay60sec/`
  - Sample `13_DMSO2_CD3OD_5mm_1728u2000uL_yuzo190815_qH_D60s_NS16-1-1.jdf`
  - EC `12_DMSO2_DMSO_5mm_1599u1500uL_yuzo190815_qH_D60s_NS16-1-1.jdf`
- D60 IC/EC calculation reproduction, FID06 rows 4–5.
  `FID06_Table4/01_Isoniazid/`
  - Sample `BTMSBandINH_5mm_1064u8463u1500uL_yuzo190817_qH_D60s_NS16-1-1.jdf`
  - EC `13_DMSO2_CD3OD_5mm_1728u2000uL_yuzo190817_qH_D60s_NS16-1-1.jdf`

## Zenodo 17083860 · CC BY 4.0

Eight complete experiments, all 1D ¹H `zgig`, D1 30 s, NS64, TD65536, SW 21.02629441 ppm (8417.5084 Hz). No missing `fid` or `acqus`. Each path below ends in `/10` under `x_qnmr_raw_data/raw data/`.

| Batch | Sample | P1 (µs) | RG |
|---|---|---:|---:|
| 1.batch | 833-1-qnmr | 10 | 42.68 |
| 1.batch | 960-qnmr | 10 | 24.32 |
| 2.batch | 832-qnmr-MK | 15 | 120.71 |
| 2.batch | 960-qnmr-MK | 15 | 48.92 |
| 3.batch | 835-qnmr | 15 | 196.86 |
| 3.batch | 961-qnmr | 15 | 39.19 |
| 4.batch | 835-2-JM-qnmr | 15 | 196.86 |
| 4.batch | 961-2-JM-qnmr | 15 | 54.82 |

**Eligibility.** Primary Bruker reader fixture. `parm.txt`, processed spectra and peak lists provide acquisition and axis cross-checks. Two probe IDs and several gain settings help test metadata checks. D1 is constant and P1 alone cannot establish flip angle. Internal and external accuracy lack a verified assay and independent reference. `intrng` and peak lists give candidate regions. Reliability labels remain open.

## Zenodo 7221753 · CC BY 4.0

715 complete 1D ¹H experiments under `x_nmr_data/NMR Data/`. All have TD65536, SW ≈20.01891529 ppm (≈8012.8205 Hz), P1 12.4 µs and D1 16 s except `quinine/9–10` at 1 s. Programs are 714 `zg` and one `zg30` (`quinine/9`). Compact NS/RG map verified against all 715 `acqus` files.

| Folder | Complete IDs | NS/RG |
|---|---|---|
| `quinine` (68) | 9,10,29,79,80 | (2,2.63), (1,206.8), (64,81.86), (512,81.86), (1024,81.86), respectively |
| `quinine` | 30–36, 40–46, 50–56, 60–66, 70–76 | Each block NS 1,2,4,8,16,32,64, RG81.86 |
| `quinine` | 131–137, 141–147, 151–157, 161–167 | Each block NS16, RG 74.8,66.74,59.25,53.19,46.47,40.22,35.61 |
| `quinine-sNN-1/-2` | NN=06 through 202, 07/08 through 171, 09 through 129, 10 through 123, all starting at 100. NN=11 has 100–105,148–153,196–202, plus 203 for `-2`. | For 100–201, successive six-ID blocks cycle NS 2,4,8,16,32,64. Block RGs from 100 are 25.35,28.45,32.29,35.61,40.22,46.47,53.19,59.25,66.74,74.8,81.86,91.38,102.91,116.93,132.74,149.03,163.23. |
| `quinine-sNN-1/-2` exceptions | 29 in all `-1`, 120 in `-1` for 06–10, 202 and 203 where present | NS/RG respectively 64/81.86, 4/35.61, 128/163.23, 128/25.35 |

No raw folder lacks `acqus`. Both `quinine-s11-1` and `-2` lack raw signal in 106–147 and 154–195 (168 folders total). Their sheets of acquisition parameters report `zg`, D1 16, NS1, TD65536, standard SW/P1 above and RG81.86. `quinine/20` and `/21` use `ser` and `acqu2s` for pseudo-2D ¹H `t1ir`, TD16384, same SW/P1, RG81.86, D1/NS 10/8 and 30/32. All seven `/29` FIDs match `quinine/36` byte for byte.

**Eligibility.** Bruker reader and diagnostic fixture, including missing FID and unsupported `t1ir`. NS/RG sweeps, processed spectra and separate inversion-recovery data support acquisition and scaling checks. Internal and external accuracy have no qualified standard for either plan model. Integral exports support reproduction. Independent region labels remain open. Apple Numbers sheets skipped.

## Case manifest

`data/case_manifest.json` records runs, preparations, dependencies, case splits, exposures, pilot pairs and frozen cohorts. All current cases are development cases.

`python tools/build_case_manifest.py --data-root PATH --out PATH` builds the snapshot from the inventory, frozen windows and local acquisition headers. Choose an output outside the data root.

`python tools/check_case_manifest.py` checks the snapshot against the committed inventory offline. It checks raw digests and processed keys. A processed digest can be checked when a run records `processed_sha256`. Keep development and final holdout cases in separate components. Excluded cases retain their links. Holdout components need known preparations and resolved exposure records, and exposed runs stay outside them.

## Redistribution status

[`data/redistribution.json`](../data/redistribution.json) lists every shipped path with its licence, source, changes and a distributable flag. `LICENSE`, data and fixture entries pin a SHA-256. Every file under `tests/fixtures` needs its own path rule. `git ls-files -z | python tools/export_check.py -z` is a dry run of an export. It exits nonzero for an unknown, missing, symlinked or changed file. It also fails when a cleared CC BY file ships without `README.md` and `tests/fixtures/README.md`, which carry its credit.

- Project code, tests, tools, docs and project made JSON fixtures: GPL-3.0-or-later. Fixtures that copy Nishizaki sheet values also carry CC0-1.0. `processing_references.json` and `integral_cases.json` also carry CC-BY-4.0 for values derived from the Zenodo 17083860 experiment.
- `LICENSE`: FSF licence text, verbatim copying only.
- Nishizaki JDF fixtures: CC0-1.0, distributable.
- Zenodo 17083860 `acqus` and `fid`: CC-BY-4.0, distributable with the credit in the root README.
- Zenodo 17083860 `acqu` and `pulseprogram`: not distributable. The reader does not need them. The pulse program is vendor text whose rights are not shown.

Open publication blockers:

- Under decision 15, build the public candidate from the cleared files with their required credits and the licence. Verify a single initial commit with Erenext LLC and the decision 12 contact for both author and committer, UTC dates and no commit trailers. Confirm that the Zenodo 17083860 `acqu` and `pulseprogram` files are absent. [docs/decisions.md](decisions.md)
- Repeat the privacy scan on the final public candidate tree and commit message. Under decision 16, keep the public-deposit fixtures byte-identical, including their original header fields. Recheck fixture checksums and required credits. [docs/open-gates.md](open-gates.md), [docs/decisions.md](decisions.md)

## Open gaps

- Qualify the JEOL pilot with the preparation records, CRM lot values and uncertainty coverage. Resolve final solution volume versus tube fill, vendor scaling, applicable T1, pulse calibration and residual transfer uncertainty.
- Keep linked stocks, repeated runs and shared EC batches together. Resolve the three weighing sessions before counting independent preparations or freezing the holdout.
- Verify axis direction and known peak positions for both vendors. Supply independently assigned regions, a Bruker same-tube accuracy set, and hand-worked and synthetic fixtures for the remaining §6 checks.
