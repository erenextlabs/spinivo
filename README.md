# Spinivo

Spinivo is a free, open source quality tool for quantitative NMR (qNMR). It reads supported raw JEOL and Bruker FIDs, checks acquisition settings and integral regions, and calculates purity with a GUM uncertainty budget when the required evidence is available. When required evidence is missing, Spinivo reports that checked purity cannot be assessed and names what is needed.

## Install and use

From a checkout of this repository, install Spinivo in a Python 3.11 or newer environment.

```sh
python -m pip install .
```

Use `inspect` for a JEOL JDF file or Bruker experiment folder. Use `assess` for a JSON assay record. Use `compare` to compare a checked assay with your own certified reference material (CRM). Each command defaults to text and accepts `--format json`. See the [assay record format](docs/assay-record.md) and the [CRM comparison](docs/crm-comparison.md).

```sh
spinivo inspect /path/to/sample.jdf --format json
spinivo assess /path/to/assay.json
spinivo compare /path/to/comparison.json
```

You can also use `python -m spinivo` in place of `spinivo`.

## Limitations

Spinivo 1.1.0 supports conventional 1D proton qNMR. A checked pass means the documented checks passed for the supplied assay. Spinivo reproduces specified calculations and GUM propagation for supported inputs. Spinivo makes no general purity accuracy claim. Budget completeness and uncertainty coverage on real samples remain unverified. The fraction of eligible preparations that yield checked purity is also unmeasured for both methods.

Reader and acquisition checks use fixtures. Integral warnings use pilot thresholds, with release scores measured on a frozen synthetic set. The baseline warning has yet to be validated independently. Spinivo reports the effect of baseline shape on the area ratio as a diagnostic and excludes it from the checked-purity budget. Nishizaki spectra and workbooks provide development evidence for both methods. Some pilot area ratios differ from workbook values.

A checked result requires signal bounds and sourced assignments, applicable relaxation and pulse calibration records, a noise region, preparation details and supported uncertainty inputs. External calibration also requires a documented pair and evidence for receiver scaling and response transfer. Spinivo withholds checked purity when required evidence is missing. Entered areas give a calculation-only value.

The `compare` command compares checked purity with the laboratory's own certified reference material (CRM) and reports signed error, uncertainty agreement and a separate verdict against the laboratory's predeclared limit ([CRM comparison](docs/crm-comparison.md)). Neither verdict qualifies a method or releases a batch. Each laboratory is responsible for its own method validation and batch release.

Release qualification uses software checks. The optional holdout study covers independent accuracy, fraction assessed and uncertainty agreement. [Software release gates](docs/open-gates.md#software-release-gates) · [Optional holdout study](docs/open-gates.md#optional-holdout-study)

## Licence

The code, tests, tools and documentation use the GNU General Public License, version 3 or later. [LICENSE](LICENSE) holds the FSF licence text, which may be copied only verbatim. Test fixtures from public deposits keep their own licences, listed below and in [`data/redistribution.json`](data/redistribution.json). Run `git ls-files -z | python tools/export_check.py -z` to list the files an export may ship.

## Data credits

- **Bruker fixture.** Jingyi Hu, "qNMRrawdata_EnhancingPolyketideProduction", Zenodo, [doi:10.5281/zenodo.17083860](https://doi.org/10.5281/zenodo.17083860), licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Spinivo ships a subset of experiment `833-1-qnmr/10`: the files `acqus` and `fid`, with unchanged bytes. Only the folder name changed. The `acqu` and `pulseprogram` files are not distributable. `tests/fixtures/processing_references.json` and `tests/fixtures/integral_cases.json` also hold values derived from this experiment.
- **JEOL fixtures and workbook values.** Nishizaki et al., Harvard Dataverse, [doi:10.7910/DVN/QFY6L9](https://doi.org/10.7910/DVN/QFY6L9), version 1.0 (depositor Guido Pauli), dedicated to the public domain under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). Two JDF files are shipped with unchanged bytes and new file names.

Some JSON fixtures include credited experimental values from the deposits above. Example uncertainties and uncertainty budgets are illustrative. [tests/fixtures/README.md](tests/fixtures/README.md)
