# People's preferences about their regional renewable electricity supply

Measuring and analysing people's preferences about renewable electricity supply within their regions.

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need [conda](https://conda.io/docs/index.html) to run the analysis. Using conda, you can create a conda environment from within you can run it:

    conda env create -f environment.yaml

## Run the analysis

    snakemake  --use-conda

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf



## Be notified of build successes or fails

  As the execution of this workflow may take a while, you can be notified whenever the execution terminates either successfully or unsuccessfully. Notifications are sent by email. To activate notifications, add the email address of the recipient to the configuration key `email`. You can add the key to your configuration file, or you can run the workflow the following way to receive notifications:

      snakemake --use-conda --config email=<your-email>

## Run the tests

    snakemake test --use-conda

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `scripts`: contains the Python source code as scripts
* `rules`: contains Snakemake rule definitions
* `envs`: contains execution environments
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`. This excludes the KlinicSlab font family (all files in `./layout/fonts/`) which is copyright Lost Type.
