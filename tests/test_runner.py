import sys

import pytest
import pandas as pd
import arviz as az


def run_test(snakemake):
    exit_code = pytest.main(
        [
            snakemake.input.test_dir,
            f"--html={snakemake.output[0]}",
            "--self-contained-html",
            "--verbose",
        ],
        plugins=[
            _create_config_plugin(snakemake)
        ]
    )
    sys.exit(exit_code)


def _create_config_plugin(snakemake):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="session")
        def config(self):
            return dict(snakemake.params.config)

        @pytest.fixture(params=snakemake.input.national_conjoints)
        def national_conjoint(self, request):
            return pd.read_feather(request.param)

        @pytest.fixture
        def conjoint(self):
            return pd.read_feather(snakemake.input.conjoint)

        @pytest.fixture
        def respondents(self, conjoint):
            return conjoint.groupby("RESPONDENT_ID").first()

        @pytest.fixture
        def country_id(self, national_conjoint):
            return national_conjoint.loc[:, "RESPONDENT_COUNTRY"].iloc[0]

        @pytest.fixture(scope="session")
        def covariate_model(self):
            return az.from_netcdf(snakemake.input.covariate_model)

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(snakemake)
