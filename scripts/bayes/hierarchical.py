from pathlib import Path

import pandas as pd
import pymc as pm
import arviz as az
import xarray as xr


ATTRIBUTES = [
    "TECHNOLOGY",
    "LAND",
    "TRANSMISSION",
    "SHARE_IMPORTS",
    "PRICES",
    "OWNERSHIP"
]
OPTIONS_PER_RESPONDENT = 16
YEAR_OF_SURVEY = 2022


def sample_prior(path_to_data: str, limit_respondents: bool, n_respondents_per_country: int,
                 model_variety: str, covariances: bool, n_draws: int, random_seed: int,
                 path_to_output: str):
    model = HierarchicalModel.for_variety(model_variety)(
        path_to_data=path_to_data,
        limit_respondents=limit_respondents,
        n_respondents_per_country=n_respondents_per_country,
        covariances=covariances
    )
    (
        pm
        .sample_prior_predictive(samples=n_draws, model=model, random_seed=random_seed)
        .to_netcdf(path_to_output)
    )


def sample_posterior(path_to_data: str, limit_respondents: bool, n_respondents_per_country: int,
                     model_variety: str, covariances: bool, n_draws: int, n_tune: int,
                     n_cores: int, random_seed: int, path_to_output: str):
    model = HierarchicalModel.for_variety(model_variety)(
        path_to_data=path_to_data,
        limit_respondents=limit_respondents,
        n_respondents_per_country=n_respondents_per_country,
        covariances=covariances
    )

    inference_data = pm.sample(
        model=model,
        draws=n_draws,
        tune=n_tune,
        cores=n_cores,
        random_seed=random_seed,
        return_inferencedata=True,
        target_accept=0.9
    )
    # The following line should not be necessary as Snakemake takes care of ensuring folders
    # exist. However, I've seen this failing on very long running jobs on the cluster in which
    # Snakemake crashed during job execution. After Snakemake's crash, the folder did not exist
    # (anymore?) and the job failed. As posterior sampling takes very long, this must not happen
    # and the following line ensures it does not.
    Path(path_to_output).parent.mkdir(exist_ok=True, parents=True)
    inference_data.to_netcdf(path_to_output)


def predict(path_to_in_sample_data: str, path_to_trace_data: str, path_to_out_sample_data: str,
            limit_respondents: bool, n_respondents_per_country: int, model_variety: str,
            covariances: bool, random_seed: int, path_to_output: str):
    if model_variety in ("mrp", "covariates"):
        raise NotImplementedError("Prediction for models with covariates is not implemented.")
    model = HierarchicalModel.for_variety(model_variety)(
        path_to_data=path_to_in_sample_data,
        limit_respondents=limit_respondents,
        n_respondents_per_country=n_respondents_per_country,
        covariances=covariances
    )
    inference_data = az.from_netcdf(path_to_trace_data)
    new_data = pd.read_feather(path_to_out_sample_data).set_index(["RESPONDENT_ID", "CHOICE_SET", "LABEL"])
    dummies = pd.get_dummies(new_data.loc[:, ATTRIBUTES], drop_first=True, prefix_sep=":")
    with model:
        pm.set_data({
            "r": new_data.xs("Left", level="LABEL").index.get_level_values("RESPONDENT_ID").remove_unused_categories().codes,
            "dummies_left": dummies.xs("Left", level="LABEL").values,
            "dummies_right": dummies.xs("Right", level="LABEL").values,
            "c": new_data.groupby("RESPONDENT_ID").RESPONDENT_COUNTRY.first().cat.codes
        })
        inference_data = pm.sample_posterior_predictive(
            inference_data,
            var_names=["p_left"],
            return_inferencedata=True,
            predictions=True,
            extend_inferencedata=False,
            random_seed=random_seed,
        )
    inference_data.to_netcdf(path_to_output)


def poststratify(path_to_conjoint: str, limit_respondents: bool, n_respondents_per_country: int,
                 model_variety: str, covariances: bool, path_to_inference_data: str, path_to_census: str,
                 path_to_output: str):
    model = HierarchicalModel.for_variety(model_variety)(
        path_to_data=path_to_conjoint,
        limit_respondents=limit_respondents,
        n_respondents_per_country=n_respondents_per_country,
        covariances=covariances
    )
    (
        model
        .poststratify(
            az.from_netcdf(path_to_inference_data),
            xr.open_dataset(path_to_census))
        .to_netcdf(path_to_output)
    )


class HierarchicalModel(pm.Model):
    variety = None
    variety_versions = {}
    covariate_col_names = []

    def __init__(self, path_to_data: str, limit_respondents: bool, n_respondents_per_country: int,
                 covariances: bool, name: str = ""):
        super().__init__(name)
        self.covariances = covariances

        self.conjoint = (
            pd
            .read_feather(path_to_data)
            .pipe(filter_respondents, limit_respondents, n_respondents_per_country)
            .pipe(self.preprocess_data) # hook for subclasses
            .apply(self.remove_unused_categories, axis="index")
            .set_index(["RESPONDENT_ID", "CHOICE_SET", "LABEL"])
        )
        self.use_imputed_covariates()
        self.respondents = self.conjoint.groupby("RESPONDENT_ID").first()
        assert self.respondents[self.covariate_col_names].notna().all(axis=None)

        dummies = pd.get_dummies(self.conjoint.loc[:, ATTRIBUTES], drop_first=True, prefix_sep=":")

        self.add_coords({
            "level": dummies.columns.values,
            "level_repeat": dummies.columns.values,
            "respondent": self.conjoint.index.get_level_values("RESPONDENT_ID").categories,
            "country": self.conjoint.RESPONDENT_COUNTRY.cat.categories,
        })

        r = pm.MutableData(
            "r",
            self.conjoint.xs("Left", level="LABEL").index.get_level_values("RESPONDENT_ID").codes,
            dims="choice_situation"
        )
        dummies_left = pm.MutableData(
            "dummies_left",
            dummies.xs("Left", level="LABEL").values,
            dims=["choice_situation", "level"]
        )
        dummies_right = pm.MutableData(
            "dummies_right",
            dummies.xs("Right", level="LABEL").values,
            dims=["choice_situation", "level"]
        )
        choice_left = pm.ConstantData(
            "choice_left",
            self.conjoint.xs("Left", level="LABEL").CHOICE_INDICATOR.values,
            dims="choice_situation"
        )
        pm.MutableData(
            "c",
            self.respondents.RESPONDENT_COUNTRY.cat.codes,
            dims=["respondent"]
        )

        partworths = self.build_partworths()

        mu_left_intercept = pm.Normal('mu_left_intercept', 0, sigma=0.25)
        sigma_left_intercept = pm.Exponential('sigma_left_intercept', 3)
        z_left_intercept = pm.Normal("z_left_intercept", 0, sigma=1, dims="respondent")
        left_intercept = pm.Deterministic( # TODO let covar with partsworths
            'left_intercept',
            mu_left_intercept + z_left_intercept * sigma_left_intercept,
            dims="respondent"
        )

        u_left = pm.Deterministic(
            "u_left",
            left_intercept[r] + pm.math.sum(partworths[r, :] * dummies_left, axis=1),
            dims="choice_situation"
        )
        u_right = pm.Deterministic(
            "u_right",
            pm.math.sum(partworths[r, :] * dummies_right, axis=1),
            dims="choice_situation"
        )
        p_left = pm.Deterministic(
            "p_left",
            pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right)),
            dims="choice_situation"
        )

        pm.Bernoulli(
            "choice",
            p=p_left,
            observed=choice_left,
            dims="choice_situation"
        )

    def preprocess_data(self, df: pd.DataFrame) -> pd.DataFrame:
        return df # generally no more preprocessing is necessary

    def remove_unused_categories(self, col: pd.Series):
        if pd.api.types.is_categorical_dtype(col):
            return col.cat.remove_unused_categories()
        else:
            return col

    def use_imputed_covariates(self):
        existing_imputed_covariates = {
            covariate: self.conjoint[f"{covariate}_imputed"]
            for covariate in self.covariate_col_names
            if f"{covariate}_imputed" in self.conjoint.columns
        }
        print(f"Using imputed data for the following covariates: {list(existing_imputed_covariates.keys())}")
        self.conjoint = self.conjoint.assign(**existing_imputed_covariates)

    def build_partworths(self):
        raise NotImplementedError("Partworth is defined for subclasses only.")

    def add_varying_effect(self, dim: str, eta: float = 4, sd: float = 2):
        z = pm.Normal(f"z_{dim}", mu=0.0, sigma=1.0, dims=["level", dim])
        if self.covariances:
            chol, rho, sigma = pm.LKJCholeskyCov(
                f"chol_{dim}",
                n=len(self.coords["level"]),
                eta=eta,
                sd_dist=pm.Exponential.dist(sd),
                compute_corr=True,
                store_in_trace=False
            )
            pm.Deterministic(f"sigma_{dim}", sigma, dims="level")
            pm.Deterministic(f"rho_{dim}", rho, dims=["level", "level_repeat"])

            effect = pm.Deterministic(f"effect_{dim}", pm.math.dot(chol, z), dims=["level", dim])
        else:
            sigma = pm.Exponential(f"sigma_{dim}", sd, dims="level")
            effect = pm.Deterministic(
                f"effect_{dim}",
                (z.T * sigma).T,
                dims=["level", dim]
            )
        return effect

    def poststratify(self, inference_data: az.InferenceData, census: xr.Dataset) -> az.InferenceData:
        raise NotImplementedError("Poststratification not implemented.")

    def poststratify_effect(self, effect: xr.DataArray, frequencies: xr.DataArray, dimension: str) -> xr.DataArray:
        """Weigh effect with frequencies of occurency along defined dimension.

        For example, takes an effect of age groups and reweighs them with the number of
        people in these age groups.
        """
        return (effect * frequencies).sum(dimension) / frequencies.sum(dimension)

    def poststratify_varying_effect(self, posterior: xr.Dataset, census: xr.Dataset, feature: str) -> xr.DataArray:
        effect = posterior[f"effect_{feature}"]
        frequencies_feature = census[f"frequency_{feature}"]
        return self.poststratify_effect(effect, frequencies_feature, feature)

    @classmethod
    def register(cls):
        HierarchicalModel.variety_versions[cls.variety] = cls

    @classmethod
    def for_variety(cls, variety_id):
        if variety_id not in HierarchicalModel.variety_versions.keys():
            raise ValueError(f"Unknown model variety: {variety_id}.")
        return HierarchicalModel.variety_versions[variety_id]


class NocovariatesModel(HierarchicalModel):
    variety = 'nocovariates'
    covariate_col_names = []

    def __init__(self, path_to_data: str, limit_respondents: bool, n_respondents_per_country: int,
                 covariances: bool, name: str = ""):
        super().__init__(path_to_data, limit_respondents, n_respondents_per_country, covariances)

    def build_partworths(self):
        alpha = pm.Normal('alpha', 0, sigma=1, dims="level")
        country = self.add_varying_effect("country", eta=4, sd=2)
        respondent = self.add_varying_effect("respondent", eta=4, sd=2)
        return pm.Deterministic(
            "partworths",
            alpha + country[:, self.c].T + respondent.T,
            dims=["respondent", "level"]
        )

    def poststratify(self, inference_data: az.InferenceData, census: xr.Dataset) -> az.InferenceData:
        posterior = inference_data.posterior
        alpha = posterior.alpha
        country_partworth = self.poststratify_effect(
            posterior.effect_country,
            census["frequency_countries"],
            "country"
        )

        partworth = (
            alpha
            + country_partworth
        )
        return az.convert_to_inference_data(
            partworth.rename("pop_means"),
            group="poststratify"
        )


class NoCovariatesVaryingVariationModel(NocovariatesModel):
    variety = 'nocovariates-varying-variation'
    covariate_col_names = []

    def __init__(self, path_to_data: str, limit_respondents: bool, n_respondents_per_country: int,
                 covariances: bool, name: str = ""):
        super().__init__(path_to_data, limit_respondents, n_respondents_per_country, covariances)

    def build_partworths(self):
        alpha = pm.Normal('alpha', 0, sigma=1, dims="level")
        country = self.add_varying_effect("country", eta=4, sd=2)

        z_respondent = pm.Normal("z_respondent", mu=0.0, sigma=1.0, dims=["level", "respondent"])
        sigma_respondent = pm.Exponential("sigma_respondent", 2, dims=["level", "country"])
        respondent = pm.Deterministic(
            "effect_respondent",
            (z_respondent.T * sigma_respondent[:, self.c].T).T,
            dims=["level", "respondent"]
        )

        return pm.Deterministic(
            "partworths",
            alpha + country[:, self.c].T + respondent.T,
            dims=["respondent", "level"]
        )


class CovariatesModel(HierarchicalModel):
    variety = 'covariates'
    covariate_col_names = [
        "Q3_GENDER",
        "Q4_BIRTH_YEAR_aggregated",
        "Q6_AREA",
        "Q9_EDUCATION"
    ]

    def __init__(self, path_to_data: str, limit_respondents: bool, n_respondents_per_country: int,
                 covariances: bool, name: str = ""):
        super().__init__(path_to_data, limit_respondents, n_respondents_per_country, covariances)

    def build_partworths(self):
        covariate_coords = {
            "gender": self.conjoint.Q3_GENDER.cat.categories,
            "age": self.conjoint.Q4_BIRTH_YEAR_aggregated.cat.categories,
            "area": self.conjoint.Q6_AREA.cat.categories,
            "education": self.conjoint.Q9_EDUCATION.cat.categories
        }
        self.add_coords(covariate_coords)

        g = pm.MutableData(
            "g",
            self.respondents.Q3_GENDER.cat.codes,
            dims="respondent"
        )
        a = pm.MutableData(
            "a",
            self.respondents.Q4_BIRTH_YEAR_aggregated.cat.codes,
            dims="respondent"
        )
        ar = pm.MutableData(
            "ar",
            self.respondents.Q6_AREA.cat.codes,
            dims="respondent"
        )
        e = pm.MutableData(
            "e",
            self.respondents.Q9_EDUCATION.cat.codes,
            dims="respondent"
        )

        alpha = pm.Normal('alpha', 0, sigma=1, dims="level")
        country = self.add_varying_effect("country", eta=4, sd=2)
        gender = self.add_varying_effect("gender", eta=4, sd=2)
        age = self.add_varying_effect("age", eta=4, sd=2)
        area = self.add_varying_effect("area", eta=4, sd=2)
        edu = self.add_varying_effect("education", eta=4, sd=2)
        respondent = self.add_varying_effect("respondent", eta=4, sd=2)
        return pm.Deterministic(
            "partworths",
            alpha + country[:, self.c].T + gender[:, g].T + age[:, a].T + area[:, ar].T + edu[:, e].T + respondent.T,
            dims=["respondent", "level"]
        )


class MrPModelAdmin0(HierarchicalModel):
    variety = 'mrp0'
    column_mapping = {
        "gender": "Q3_GENDER",
        "age": "Q4_BIRTH_YEAR_aggregated",
        "education": "Q9_EDUCATION",
        "concern": "Q11_CLIMATE_CONCERN",
        "justice": "Q15_ATTRIBUTE_IMPORTANCE_O1",
        "effectiveness": "Q22_CLIMATE_PROTECTION"
    }
    index_mapping = {
        "gender": "g",
        "age": "a",
        "education": "e",
        "concern": "cc",
        "justice": "j",
        "effectiveness": "eff"
    }
    covariate_col_names = column_mapping.values()
    covariates = column_mapping.keys()

    def __init__(self, path_to_data: str, limit_respondents: bool, n_respondents_per_country: int,
                 covariances: bool, name: str = ""):
        assert set(self.covariate_col_names) == set(self.column_mapping.values())
        assert set(self.column_mapping.keys()) == set(self.index_mapping.keys())
        super().__init__(path_to_data, limit_respondents, n_respondents_per_country, covariances)

    def preprocess_data(self, df: pd.DataFrame) -> pd.DataFrame:
        # Remove respondents without location information if region is a covariate.
        # Removes 48 (1.2%) people in total.
        # Removes 3 (0.3%) Germans.
        # Removes 4 (0.4%) Danes.
        # Removes 29 (2.8%) Poles.
        # Removes 12 (1.2%) Portuguese.
        if "region" in self.covariates:
            df = (
                df
                .dropna(axis="index", subset=[self.column_mapping["region"]])
            )
        return df

    def build_partworths(self):
        covariate_coords = {
            covariate_name: self.conjoint[column_name].cat.categories
            for covariate_name, column_name in self.column_mapping.items()
        }
        self.add_coords(covariate_coords)

        indices = [
            pm.MutableData(
                index,
                self.respondents[self.column_mapping[covariate]].cat.codes,
                dims="respondent"
            )
            for covariate, index in self.index_mapping.items()
        ]

        alpha = pm.Normal('alpha', 0, sigma=1, dims="level")
        country = self.add_varying_effect("country", eta=4, sd=2)
        respondent = self.add_varying_effect("respondent", eta=4, sd=2)

        effects = [
            self.add_varying_effect(covariate, eta=4, sd=2)
            for covariate in self.index_mapping.keys()
        ]

        return pm.Deterministic(
            "partworths",
            (
                alpha + country[:, self.c].T + respondent.T
                + sum(effect[:, index].T for effect, index in zip(effects, indices))
            ),
            dims=["respondent", "level"]
        )

    def poststratify(self, inference_data: az.InferenceData, census: xr.Dataset) -> az.InferenceData:
        national = self.poststratify_national(inference_data, census)
        sample = self.sample_partworths(inference_data, census)
        bias = national - sample

        all = xr.Dataset({
            "national_partworths": national,
            "sample_partworths": sample,
            "national_bias": bias
        })
        return az.convert_to_inference_data(
            all,
            group="poststratify"
        )

    def poststratify_national(self, inference_data: az.InferenceData, census: xr.Dataset) -> az.InferenceData:
        posterior = inference_data.posterior
        national_census = (
            census
            .drop("frequency_countries") # TODO this should not be there in the first place
            .drop_dims("country") # TODO this should not be there in the first place
            .groupby("country_of_admin1")
            .sum()
            .rename(country_of_admin1="country")
        )

        alpha = posterior.alpha
        country_partworth = posterior.effect_country

        covariate_partworths = [
            self.poststratify_varying_effect(posterior, national_census, covariate)
            for covariate in set(self.covariates) - set(["region"])
        ]

        partworth = (
            alpha
            + country_partworth
            + sum(covariate_partworths)
        )
        return partworth

    def sample_partworths(self, inference_data: az.InferenceData, census: xr.Dataset) -> az.InferenceData:
        posterior = inference_data.posterior
        constant_data = inference_data.constant_data

        alpha = posterior.alpha
        country_partworth = self.sample_based_feature(posterior, "country", constant_data.c)

        covariate_partworths = [
            self.sample_based_feature(posterior, covariate, constant_data[index])
            for covariate, index in self.index_mapping.items()
        ]

        partworth = (
            alpha
            + country_partworth
            + sum(covariate_partworths)
        )
        return partworth

    def sample_based_feature(self, posterior: xr.Dataset, feature: str,
                             respondent_features: xr.DataArray) -> xr.DataArray:
        """Determine effect size based on sample characteristics rather than poststratification."""
        effect = posterior[f"effect_{feature}"]
        return effect.isel({feature: respondent_features}).mean("respondent")


class MrPModelAdmin1(MrPModelAdmin0):
    variety = 'mrp1'
    column_mapping = {
        "region": "RESPONDENT_ADM1"
    }
    index_mapping = {
        "region": "rgn"
    }
    covariate_col_names = column_mapping.values()
    covariates = column_mapping.keys()

    def poststratify(self, inference_data: az.InferenceData, census: xr.Dataset) -> az.InferenceData:
        national = self.poststratify_national(inference_data, census)
        regional = self.poststratify_regional(inference_data, census)
        sample = self.sample_partworths(inference_data, census)
        bias = national - sample

        all = xr.Dataset({
            "national_partworths": national,
            "regional_partworths": regional,
            "sample_partworths": sample,
            "national_bias": bias
        })
        return az.convert_to_inference_data(
            all,
            group="poststratify"
        )

    def poststratify_regional(self, inference_data: az.InferenceData, census: xr.Dataset) -> az.InferenceData:
        posterior = inference_data.posterior
        alpha = posterior.alpha
        country_partworth = posterior.effect_country
        region_partworth = posterior["effect_region"]

        covariate_partworths = [
            self.poststratify_varying_effect(posterior, census, covariate)
            for covariate in set(self.covariates) - set(["region"])
        ]

        partworth = (
            alpha
            + country_partworth
            + region_partworth
            + sum(covariate_partworths)
        )
        return partworth


class MrPModelAdmin2(MrPModelAdmin1):
    variety = 'mrp2'
    column_mapping = {
        "region": "RESPONDENT_ADM2"
    }
    index_mapping = {
        "region": "rgn"
    }
    covariate_col_names = column_mapping.values()
    covariates = column_mapping.keys()


NocovariatesModel.register()
NoCovariatesVaryingVariationModel.register()
CovariatesModel.register()
MrPModelAdmin0.register()
MrPModelAdmin1.register()
MrPModelAdmin2.register()


def filter_respondents(df, limit_respondents, n_respondents_per_country):
    if limit_respondents:
        df = df.groupby("RESPONDENT_COUNTRY").head(n_respondents_per_country * OPTIONS_PER_RESPONDENT)
    return df


if __name__ == "__main__":
    match snakemake.wildcards.sample:
        case "prior":
            sample_prior(
                path_to_data=snakemake.input.data,
                n_draws=snakemake.params.n_draws,
                limit_respondents=bool(snakemake.params.limit_respondents),
                n_respondents_per_country=int(snakemake.params.limit_respondents) // 4,
                random_seed=snakemake.params.random_seed,
                model_variety=snakemake.params.model_variety,
                covariances=snakemake.params.covariances,
                path_to_output=snakemake.output[0]
            )
        case "posterior":
            sample_posterior(
                path_to_data=snakemake.input.data,
                n_tune=snakemake.params.n_tune,
                n_draws=snakemake.params.n_draws,
                limit_respondents=bool(snakemake.params.limit_respondents),
                n_respondents_per_country=int(snakemake.params.limit_respondents) // 4,
                n_cores=snakemake.threads,
                random_seed=snakemake.params.random_seed,
                model_variety=snakemake.params.model_variety,
                covariances=snakemake.params.covariances,
                path_to_output=snakemake.output[0]
            )
        case "prediction":
            predict(
                path_to_in_sample_data=snakemake.input.in_sample,
                path_to_out_sample_data=snakemake.input.out_sample,
                path_to_trace_data=snakemake.input.trace,
                limit_respondents=bool(snakemake.params.limit_respondents),
                n_respondents_per_country=int(snakemake.params.limit_respondents) // 4,
                random_seed=snakemake.params.random_seed,
                model_variety=snakemake.params.model_variety,
                covariances=snakemake.params.covariances,
                path_to_output=snakemake.output[0]
            )
        case "poststratify":
            poststratify(
                path_to_conjoint=snakemake.input.conjoint,
                limit_respondents=snakemake.params.limit_respondents,
                n_respondents_per_country=int(snakemake.params.limit_respondents) // 4,
                model_variety=snakemake.params.model_variety,
                covariances=snakemake.params.covariances,
                path_to_inference_data=snakemake.input.inference_data,
                path_to_census=snakemake.input.census,
                path_to_output=snakemake.output[0]
            )
        case _:
            raise ValueError(f"Unknown type {snakemake.wildcards.sample}.")
