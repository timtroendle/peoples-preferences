import sys

import pytest
import pandas as pd


def run_test(path_to_test_dir: str, paths_to_national_conjoints: list[str], path_to_output: str, config: dict):
    exit_code = pytest.main(
        [
            path_to_test_dir,
            f"--html={path_to_output}",
            "--self-contained-html",
            "--verbose",
        ],
        plugins=[
            _create_config_plugin(
                config=config,
                paths_to_national_conjoints=paths_to_national_conjoints,
            )
        ]
    )
    sys.exit(exit_code)


def _create_config_plugin(config: dict, paths_to_national_conjoints: list[str]):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="session")
        def config(self):
            return dict(config)

        @pytest.fixture(params=paths_to_national_conjoints)
        def national_conjoint(self, request):
            return pd.read_feather(request.param)

        @pytest.fixture
        def country_id(self, national_conjoint):
            return national_conjoint.loc[:, "RESPONDENT_COUNTRY"].iloc[0]

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(
        path_to_test_dir=snakemake.input.test_dir,
        paths_to_national_conjoints=snakemake.input.national_conjoints,
        config=snakemake.params.config,
        path_to_output=snakemake.output[0]
    )
