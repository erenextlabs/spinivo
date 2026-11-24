# Test fixtures

The public release includes these unmodified raw NMR files from public deposits for reader tests. Only the file or folder names changed.

| Fixture | Source | Original path in the deposit | Licence | SHA-256 |
|---|---|---|---|---|
| `jeol/nishizaki_fid02_ma_dmso2_d20s_ns32.jdf` | Nishizaki et al., Harvard Dataverse, doi:10.7910/DVN/QFY6L9, version 1.0 (depositor Guido Pauli) | `334_Dataverse_20200705/FID02_Figure2/01_5mmNMRtube/Pattern1/0619_2019/01_MAandDMSO2_T1A_5mm_4226u1534u6000uL_yuzo190619_RG30_No90pw_D20s_NS32_qH-1-1.jdf` | CC0-1.0 | `dd9aeb9cafebaa99755ab0576bd8e1c3ac026970766db58e7a777732cd0bcaf3` |
| `jeol/nishizaki_fid01_pulse_width_array.jdf` | same deposit | `334_Dataverse_20200705/FID01_Table2/Day2/Delay120sec/MAandCAF_T1A_5mm_6435u15130u9000uL_yuzo19807_Delay120_CAF_H1-1-1.jdf` | CC0-1.0 | `fc9b1dd417bcc15912cacd655cb3891116e5baf26489c05a0e8dac94e9f483f5` |
| `bruker/zenodo17083860_833-1/10/acqus` | Jingyi Hu, "qNMRrawdata_EnhancingPolyketideProduction", Zenodo, doi:10.5281/zenodo.17083860 | `raw data/1.batch/833-1-qnmr/10/acqus` | CC-BY-4.0 | `3aed337ada7a9682461dadf1baab4ee5a3c4dbd55bf221062d285c600c1b11f2` |
| `bruker/zenodo17083860_833-1/10/fid` | same deposit | `raw data/1.batch/833-1-qnmr/10/fid` | CC-BY-4.0 | `52fae41b0e075ed323fab8088b87d6a78081f6cb4ed3bdb89cc0171e772b9251` |

Facts checked on 22 Nov 2026 with nmrglue 0.12 (independent of Spinivo code):

- JEOL 1D fixture: 1D complex FID, 65,536 points, `x_sweep` 7494.0048 Hz, `x_freq` 399782198.37825 Hz, `x_offset` 5.0 ppm, `relaxation_delay` 20.0 s, `x_pulse` 7.272 µs, `x_90_width` 7.272 µs, `x_angle` 90°, `scans` 32, `x_prescans` 4, `recvr_gain` 30, `temp_get` 25 °C, `experiment` `single_pulse_dec`, `digital_filter` TRUE. Same-tube maleic acid and dimethyl sulfone.
- JEOL array fixture: `small_two_d` layout, 1024 × 12 points, second axis `x_90_width` (pulse width array). It is not a 1D quantitative FID.
- Bruker fixture: `NUC1` 1H, `TD` 65536 (32,768 complex points), `SW` 21.0262944062838 ppm, `SW_h` 8417.50841750842 Hz, `SFO1` 400.33247204 MHz, `BF1` 400.33 MHz, `O1` 2472.04 Hz, `D1` 30 s, `P1` 10 µs, `NS` 64, `DS` 4, `RG` 42.68, `PULPROG` `zgig`, `GRPDLY` 67.9872589111328, `DECIM` 2376, `DSPFVS` 20, `TE` 298.5974 K, `NC` -6, `BYTORDA` 0, `DTYPA` 0. (Corrected 22 Nov 2026: the first version listed P1, RG and TE of another acquisition in the deposit.)

Reader tests build damaged or edited copies under a temporary directory. Bruker cases include cut files, directories with `ser` alone or alongside `fid`, conflicting dimension metadata, `D` and `P` array count mismatches or invalid tokens, and an undecodable byte in `PULPROG` after valid UTF-8 or Windows-1252 text or an intervening JCAMP comment. Both readers distinguish absent or blank nucleus metadata from an explicit non-proton value. Bruker also tests a trailing comment and the `off` value in `NUC1`. These cases use synthetic edits to the committed fixtures. There is no physical non-proton or real 2D `ser` acquisition among the fixtures.

Attribution for the Bruker files (CC BY 4.0): Jingyi Hu, qNMRrawdata_EnhancingPolyketideProduction, Zenodo, https://doi.org/10.5281/zenodo.17083860.

The Bruker values derived from this experiment in `processing_references.json` and `integral_cases.json` carry the same credit. Only `acqus` and `fid` are distributable. `acqu` and the vendor `pulseprogram` stay out of exports; the reader needs neither. Every fixture's licence and SHA-256 is in [`data/redistribution.json`](../../data/redistribution.json). Update its hash there after any fixture change.

## Integral labelling fixtures

[`integral_cases.json`](integral_cases.json) records the labels and evidence under the [pilot rubric](../../docs/integral-rubric.md). Synthetic recipes and their generator use the repository licence, GPL-3.0-or-later. The generator recreates the ideal target-plus-anchors FID, corruption and fixed-seed complex Gaussian noise. Synthetic bounds are in Hz relative to the carrier. The high-S/N recipes use nominal S/N 100,000 with baseline exclusions extending 250 Hz around the target and each phase anchor. These exclusions keep measurable Lorentzian wings out of the baseline fit while preserving the scored integration bounds.

Real bounds are in ppm. `edge_area_relative_changes` holds signed fractions in this order: lower edge inward, lower edge outward, upper edge inward, upper edge outward. Each trial moves one edge by `edge_step_ppm`. Committed files retain the licences above. The external FID05 path is relative to the extracted Nishizaki CC0-1.0 deposit (doi:10.7910/DVN/QFY6L9). Set `SPINIVO_NISHIZAKI_ROOT` to the directory containing `334_Dataverse_20200705` to run the cached-file test. The suite always checks the external path against the deposit inventory.

## Internal standard calculation fixture

The isoniazid and BTMSB assay uses Sheet1 row 5 of the Nishizaki FID06 workbook. The sheet values use CC0-1.0. The illustrative standard uncertainties, weighing correlation, coverage factor and hand-worked budget use GPL-3.0-or-later. [Fixture](internal_standard_cases.json) · [Deposit](https://doi.org/10.7910/DVN/QFY6L9)

## External calibration calculation fixture

The caffeine 04 and dimethyl sulfone 01 pair uses rows 149 and 146 of the Nishizaki FID02 workbook. The sheet values use CC0-1.0. The run and preparation IDs are fixture labels. The uncertainties, correlations, coverage factor and drift factor are illustrative and use GPL-3.0-or-later, as does the hand-worked budget. [Fixture](external_calibration_cases.json) · [Pilot pair](../../docs/pilot-pairs.md#1-fid02--caffeine-04--dmso2-01) · [Deposit](https://doi.org/10.7910/DVN/QFY6L9)
