# Spinivo

Spinivo is a free, open source quality tool for quantitative NMR (qNMR). It reads a raw JEOL or Bruker FID, checks the acquisition settings, flags unreliable integrals and computes purity with a GUM uncertainty budget. When the evidence is missing, it reports that the result cannot be assessed.

## Install and use

From a checkout of this repository, install Spinivo in a Python 3.11 or newer environment.

```sh
python -m pip install .
```

Use `inspect` for a JEOL JDF file or Bruker experiment folder. Use `assess` for a JSON assay record. Both commands default to text and accept `--format json`. See the [assay record format](docs/assay-record.md).

```sh
spinivo inspect /path/to/sample.jdf --format json
spinivo assess /path/to/assay.json
```

You can also use `python -m spinivo` in place of `spinivo`.

## Limitations

Spinivo v1.0 is planned as a limited preview for conventional 1D proton qNMR from JEOL and Bruker raw data. Both the internal-standard and external-calibration methods still need an independent test set. Pilot spectra and workbooks help check file handling and calculations. Accuracy, the fraction of preparations that yield checked purity, and agreement within reported uncertainty remain unverified for both methods.

Reader and acquisition checks have fixture coverage. Further validation is needed for both vendors. Integral warnings use pilot thresholds and need testing on independently labelled real spectra. The baseline warning needs independent validation. Spinivo reports the effect of baseline shape on the area ratio separately from the checked-purity budget. External calibration still needs a verified comparison of receiver scaling and sample response across tubes. Some pilot area ratios differ from workbook values.

A checked result requires signal bounds and sourced assignments, applicable relaxation and pulse calibration records, a noise region, preparation details and supported uncertainty inputs. External calibration also requires a documented pair and evidence for the relative receiver scale. Spinivo withholds checked purity when required evidence is missing. Entered areas give a calculation only. See the [open gates](docs/open-gates.md) for the evidence still needed and the conditions for publication.

## Licence

The code, tests, tools and documentation use the GNU General Public License, version 3 or later. [LICENSE](LICENSE) holds the FSF licence text, which may be copied only verbatim. Test fixtures from public deposits keep their own licences, listed below and in [`data/redistribution.json`](data/redistribution.json). Run `git ls-files -z | python tools/export_check.py -z` to list the files an export may ship.

## Data credits

- **Bruker fixture.** Jingyi Hu, "qNMRrawdata_EnhancingPolyketideProduction", Zenodo, [doi:10.5281/zenodo.17083860](https://doi.org/10.5281/zenodo.17083860), licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Spinivo ships a subset of experiment `833-1-qnmr/10`: the files `acqus` and `fid`, with unchanged bytes. Only the folder name changed. The `acqu` and `pulseprogram` files are not distributable. `tests/fixtures/processing_references.json` and `tests/fixtures/integral_cases.json` also hold values derived from this experiment.
- **JEOL fixtures and workbook values.** Nishizaki et al., Harvard Dataverse, [doi:10.7910/DVN/QFY6L9](https://doi.org/10.7910/DVN/QFY6L9), version 1.0 (depositor Guido Pauli), dedicated to the public domain under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). Two JDF files are shipped with unchanged bytes and new file names.

Some JSON fixtures include credited experimental values from the deposits above. Example uncertainties and uncertainty budgets are illustrative. [tests/fixtures/README.md](tests/fixtures/README.md)
