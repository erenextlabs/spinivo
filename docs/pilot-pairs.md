# JEOL pilot pairs

Deposit `10.7910/DVN/QFY6L9`. File IDs refer to the `-1-1.jdf` runs in `docs/data-inventory.md`. Paired values are sample / EC. NF means **NOT FOUND**. Concentrations are in mmol/L. Δ is recalculated minus cached Excel value in percentage points.

All six JDF headers show 1D ¹H, probe 3448, 399.782198 MHz, four prescans, 8.745124 s acquisition and nominal 90° excitation. Each has a `-1-2.jdf` counterpart with the same acquisition end time and pulse width and a different stored data array. Exact file suffixes and integration bounds are NF in the sheets.

The purity preflight checks the case manifest and permits only these six raw inventory keys. External calibration requires the documented sample and EC order. The reference comparison checks every pair before reading any pilot spectrum.

## 1. FID02 · Caffeine 04 / DMSO2 01

`FID02_CalculationSpreadsheet.xlsx`, **5-mm NMR tubes**. A124 identifies Pattern 4, A146 = `0724_2019`, and B149/B146 = 4/1.

| Quantity | Filename | Sheet | Header | Match |
|---|---|---|---|---|
| MA + solute, volume | 6.435 + 15.130 mg, 9000 µL / 4.226 + 1.534 mg, 6000 µL | C149/E149 = 6.153915741/8.648385602. C146/D146 = 6.062109072/2.713385743 | Same recipe tokens in titles | Concentrations agree |
| D1 / NS / gain | 20 s / NF / 30 each | NF | 20 s / 32 / 30 each | Yes |
| `x_pulse` (µs) | NF | F149/F146 = 7.36572/7.37538 | 7.36572/7.37538 | Sheet agrees |
| Points / SW (Hz) | NF | NF | 65536 / 7494.005 each, clipped 5995.204 | Yes |
| Decoupling / temperature | NF | NF | ¹³C MPF8 / 25 °C each | Yes |

Both runs use DMSO-d6. Direct mass and volume cells are NF. The concentrations reproduce the recipes using a 0.999 factor and inferred molar masses of 116.07 for MA, 194.19 for caffeine and 94.13 for DMSO2 in g/mol.

- Areas G149/I149 = 0.26387499/0.55648986 (MA/CAF_H7). G146/H146 = 0.26193853/0.35070252 (MA/DMSO2).
- EC N149 `=D146*(I149/H146)*(6/3)*(F149/F146)` = 8.599846658. O149 `=N149/E149` = 0.9943875139, or **99.44%** relative (Δ −7.0 × 10⁻¹⁵ pp). The pulse factor is sample/EC. Volumes enter through concentration.
- Same-tube MA K149 `=C149*(I149/G149)*(2/3)` = 8.652055806. L149 `=K149/E149` = 1.0004243803, or **100.0%** relative (Δ +4.7 × 10⁻¹⁵ pp). On the inferred 99.9% caffeine basis, EC and MA give 99.34% and 99.94% mass fraction.
- EC tube’s own MA result K146 `=C146*(H146/G146)*(2/6)` = 2.705465449. L146 `=K146/D146` = 0.9970810291.
- Alternate ECs also use the caffeine run. Q149/T149 follow `=D{r}*(I149/H{r})*(6/3)*(F149/F{r})` for r = 147/148 and give 8.651580121/8.608810928. R149 `=Q149/E149` and U149 `=T149/E149` give 1.000369378/0.9954240392.

**Gaps**

- Balance and volume uncertainties, weighing-session link, lot certificates and coverage for caffeine/DMSO2/MA, and fitted signal T1 values are NF. FID01 has a caffeine/MA delay series with the same recipe and a different header sample number.

**QUALIFIED WITH GAPS** for step 1 EC and same-tube MA calculation pilots. The run settings and workbook mapping support reproduction.

## 2. FID05 · DMSO2 13 / DMSO2 12′

`FID05_CalculationSpreadsheet.xlsx`, **Sheet1**. A5 = N1, B7/B6 = 13/12′ in the decoupled block.

| Quantity | Filename | Sheet | Header | Match |
|---|---|---|---|---|
| DMSO2, volume | 1.728 mg, 2000 µL / 1.599 mg, 1500 µL | E7/E6 = 9.169616488/11.313438861 | Same recipe tokens in titles | Concentrations agree |
| Solvent | CD3OD / DMSO | C7/C6 = CD3OD/DMSO-d6 | Methanol-d4 / DMSO-d6 | Records agree, solvents differ |
| D1 / NS / gain | 60 s / 16 / NF each | D7/D6 = 60/60 s. NS and gain NF | 60 s / 16 / 30 each | Yes |
| `x_pulse` (µs) | NF | F7/F6 = 7.27788/7.30474 | 7.27788/7.30474 | Sheet agrees |
| Points / SW (Hz) | NF | NF | 65536 / 7494.005 each, clipped 5995.204 | Yes |
| Decoupling / temperature | NF | A1 ¹³C decoupled. Temperature NF | ¹³C MPF8 / 25 °C each | Yes |

Direct mass and volume cells are NF. Both concentrations fit the recipes at 0.999 and 94.13 g/mol DMSO2.

- Areas G7/G6 = 1.19514102/1.48424385. J7 `=E6*(G7/G6)*(6/6)*(F7/F6)` = 9.076296123. K7 `=J7/E7` = 0.9898228716, or **98.98%** relative (Δ +1.1 × 10⁻¹⁴ pp). Volumes enter through concentration. On the inferred 99.9% basis this is 98.88% mass fraction.
- The sample also uses D20 EC 12. D5/E5/F5/G5 = 20/11.313438861/7.31322/1.44465037. H7 `=E5*(G7/G5)*(6/6)*(F7/F5)` = 9.314237035. I7 `=H7/E7` = 1.015771712.
- The D60 EC also feeds J8–J10 `=E6*(G{r}/G6)*(6/6)*(F{r}/F6)` = 10.32021048, 10.02572908, 12.51440927. K8–K10 `=J{r}/E{r}` = 1.017169274, 1.017959639, 1.005251922, for r = 8, 9, 10.

**Gaps**

- Balance and volume uncertainties, weighing-session link, DMSO2 lot certificate and coverage, and fitted T1 values are NF. FID05 `03_T1` has a raw T1 array for 13 (No71). A corresponding array for 12 is NF.

**NOT QUALIFIED** for matched EC purity. The solvent changes and validation of the pulse-width correction for this transfer is NF. The pair can serve as a solvent-transfer challenge.

## 3. FID06 · INH + BTMSB / DMSO2 13

`FID06_CalculationSpreadsheet.xlsx`, **Sheet1**. A5 = `INH + IC`, A4 = `13 (EC)`.

| Quantity | Filename | Sheet | Header | Match |
|---|---|---|---|---|
| Mass, volume | 1.064 mg BTMSB + 8.463 mg INH, 1500 µL / 1.728 mg DMSO2, 2000 µL | C5/D5 use exact sample amounts. B4 = 9.169616488 | Same recipe tokens in titles | Sample exact, EC consistent |
| D1 / NS / gain | 60 s / 16 / NF each | NF | 60 s / 16 / 30 each | Yes |
| `x_pulse` (µs) | NF | E5/E4 = 7.30976/7.29494 | 7.30976/7.29494 | Sheet agrees |
| Points / SW (Hz) | NF | NF | 65536 / 7494.005 each, clipped 5995.204 | Yes |
| Decoupling / temperature | NF | NF | ¹³C MPF8 / 25 °C each | Yes |

Both runs use methanol-d4. C5 `=(1.064*0.999)/226.5/1.5*1000` = 3.128582781. D5 `=8.463/137.139/1.5*1000` = 41.14074042. A10 gives 99.9% for BTMSB and DMSO2. B4 fits the EC recipe at 0.999 and 94.13 g/mol. Direct EC mass and volume cells are NF.

- Areas F4/G5/H5/I5 = 1.19608545/1.23461854/1.79906664/1.80570143.
- EC N5 `=B4*((H5+I5)/F4)*(6/4)*(E5/E4)` = 41.53736522. O5 `=N5/D5` = 1.009640682, or **101.0%** (Δ −2.3 × 10⁻¹⁴ pp). Volumes enter through concentration.
- IC K5 `=C5*((H5+I5)/G5)*(18/4)*(E5/E5)` = 41.10595076. L5 `=K5/D5` = 0.9991543747, or **99.92%** (Δ −2.9 × 10⁻¹⁴ pp). The pulse factor is unity and the shared volume cancels.

**Gaps**

- Balance and volume uncertainties, weighing-session link, lot certificates and coverage, independent INH purity and fitted T1 values are NF. EC 13 shares header sample No71 and the recipe with pair 2’s sample, with a distinct acquisition.

**QUALIFIED WITH GAPS** for step 1 IC/EC calculation reproduction. Independent INH accuracy and the GUM budget remain open.

For all pairs, tube fill, validated receiver scaling, residual response-transfer uncertainty and integral repeatability are NF. Flip-angle validation, quantitative recovery and preparation independence remain unverified.

## Stage diagnostics

Run `python tools/pilot_stages.py DATA_ROOT` with the local cache or extracted deposit directory. The tool checks the pinned pairs before reading any FID. The baseline uses the same settings and windows as `compare_references.py`. Each trial starts from that baseline. The last column is `trial ratio / baseline ratio − 1` as a percent. Same-tube MA controls appear above the cross-tube ratios.

The scan-count trial assumes stored sums. Receiver response, workbook intensity scaling, workbook ppm shift and workbook bounds are unknown. The axis trial changes FFT sampling density. Each region trial moves one window's upper edge by the larger of one acquired ppm step and 5% of the window width on the baseline trace. Pulse trials multiply the raw area ratio by the two reciprocal nominal width ratios for arithmetic comparison. Physical transfer remains unknown.

Workbook area ratios provide reproduction context. The committed evidence supports sensitivity comparisons and leaves the EC validation gate closed. FID05 is a solvent-transfer challenge. FID06 is a calculation reproduction case.
