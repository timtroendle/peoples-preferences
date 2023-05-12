# Method

## Experimental design

To assess citizen preferences for renewable electricity supply in their region, we conducted a choice experiment. Choice experiments are a widely used method to assess stated preferences in social and political science [@Hainmueller:2014; @Bansak:2021]. In our experiment, respondents are repeatedly presented with a choice between two hypothetical designs of the electricity supply in their region. Each option is formed out of six attributes describing the electricity supply: dominant technology, land requirements, level of electricity imports into the region, household electricity prices, overhead transmission capacity expansion, and ownership of the assets (see Table XX for all attributes and attribute levels; see Figure YY for an example of the choice display). Each participant was presented with a total of eight such choices from which we identify the relative importance of the attributes and their 25 levels. We randomised the combinations of attribute levels (fully randomised design, @fig:experimental-design) and the order in which attributes are presented (across respondents, but not within respondent).

The choice experiment was set up through conjointly.com and distributed with support from polling agencies in each country. Each respondent was presented with eight consecutive pairs of hypothetical regional system designs (in the following profile(s)) for a future fully renewable electricity supply with the task to select between two options. The display of attribute levels was enhanced with small pictograms for enhancing the understanding (see Figure XX).

## Sampling

We conducted our choice experiment as an online survey between 24 January and 8 February 2022 using the platform conjointly.com. The sample was drawn from opt-in panels in the four European countries using the panel service of bildendi (respondi) and partners in Germany, Poland, Denmark and Portugal. We used a hard quota on the gender (50%), and a soft quota on respondents from rural areas (>=30%).

## Country selection

FROM FRANZISKA

## Data analysis

We derive preferences from recorded choices using multinomial logit hierarchical bayes, a method commonly used for choice experiments [@Green:2004]. It is based on the random utility theory and assumes that each option has a distinct utility to each respondent and that respondents choose the option with higher utility. As we did not measure utility but choices, utility is a latent variable that is estimated by the model. The full model is shown in @eq:model.

Given that we have exactly two options per task, we model choices as a Bernoulli variable (first line in @eq:model). The probability of the Bernoulli variable depends on the deterministic parts of the utilities of both options (left option and right option, second line in @eq:model). Utilities (V) are linear combinations of the partworth utilities of each level that is included in an option (third and forth line in @eq:model). For example, if the left option showed an electricity supply based mainly on publicly-owned wind turbines within the region, utility is the sum of the partworth utilities of the attribute levels of "public utility" and "wind turbine", and the four other attribute levels included in the option. To capture the possibility that respondents preferred the left or right option irrespective of the shown attribute levels, we add an intercept term to the left utility (third line in @eq:model). The intercept turns out to be very small or zero (@fig:left-intercept) suggesting that there is no impactful bias towards left or right options. Finally, the partworth utility of each level is the sum of level-specific intercept and varying (random) effects. In the base model, we add a varying effect for country (n=4) and a varying effect for each respondent (n=4103).

$$
\begin{array}{rcl}
\text{choice}_\text{left} &\sim& \operatorname{Bern}(
p_\text{left})\\
p_\text{left} &=& \frac{\exp(V_\text{left})}{\exp(V_\text{left}) + \exp(V_\text{right})}\\
V_\text{left} &=& \alpha_\text{left} + \sum_\text{level}{x_\text{level} * \beta_\text{level}}\\
V_\text{right} &=& \sum_\text{level}{x_\text{level} * \beta_\text{level}}\\
\beta_\text{level} &=& \alpha + \beta_\text{Country} + \beta_\text{Respondent}
\end{array}
$$ {#eq:model}

Being a Bayesian model, we add prior probabilities for all parameters (@eq:priors). We use weakly informative priors to avoid unrealistically large parameter values. Given that utilities are defined on the logit scale in this model, a utility value of 4 or -4 means that an option is chosen or rejected with a probability larger than 98% (when the utility of the other is 0). Therefore, we deem partworth utilities with absolute value larger than four unrealistic and tune the prior probabilities accordingly. We model varying effects of country and respondent with mean zero (second line in @eq:priors) as the models includes a separate intercept term per partworth utility (first line in @eq:priors, bottom line in @eq:model). We do not model covariances between attribute levels as the additional computational complexity is restrictive.

$$
\begin{array}{rcl}
\alpha &\sim& \operatorname{N}(0, ~1)\\
\beta_\text{...} &\sim& \operatorname{N}(0, ~\sigma_{...})\\
\sigma_{...} &\sim& \operatorname{Exp}(2)\\
\alpha_\text{left} &\sim& \operatorname{N}(\mu_\text{left}, ~\sigma_\text{left})\\
\mu_\text{left} &\sim& \operatorname{N}(0, 0.25)\\
\sigma_\text{left} &\sim& \operatorname{Exp}(3)\\
\text{...} &\in& [\text{Country, Respondent}]\\
\end{array}
$$ {#eq:priors}

We implement our probabilistic model using PyMC [@ThomasWiecki:2023]. Code and data to reproduce our analysis are publicly available (after publication of the article). We sample the posterior distribution of all model parameters using a NUTS sampler. We run four independent chains of the Markov chain Monte Carlo to check for convergence. Each chain iterates a total of 4,000 times of which 2,000 iterations are tuning steps which we discard. The chains converge to the posterior distributions (@tbl:sample-statistics). Calculations were carried out on the ETH Euler cluster.

In addition to the base model, we implement a model that includes respondent-level covariates as varying effects. We add age (n=7), gender (n=3), education (n=7), and area (urban/rural, n=2) varying effects to better estimate the potential bias introduced through non-random sampling. In this covariate model, all additional varying effects are added to the partworth utilties of each attribute level (@eq:model2). We find that the impact of these respondent-level covariates are small and, with it, the bias introduced through non-random sampling (@fig:max-subgroup-effect). Therefore, we exclude these covariates in the base model.

$$
\beta_\text{level} = \alpha + \beta_\text{Country} + \beta_\text{Age} + \beta_\text{Gender} +\beta_\text{Education} + \beta_\text{Area} + \beta_\text{Respondent}
$$ {#eq:model2}

While the data are complete for the experimental variables, there are missing values for the covariates. We are handling these missing values for the covariate model (@eq:model2). There are missing values in the age, education, and area covariates, for which 17, 40, and 54 respondents respectively did not state valid values (0.4%, 1.0%, 1.3%). We use single, multivariate data imputation  [@Stekhoven:2012] to fill the missing values. The out-of-bag error of the imputation is 0.01, 0.32, and 0.35, respectively for the age, education, and area covariates.

# Supplementary Material

## Respondents

```table
---
caption: 'Statistics of respondents in Germany'
alignment: LRRRRRR
include: build/results/analysis/respondent-stats-DEU.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.4
    - 0.2
---
```

<div class="pagebreak">

```table
---
caption: 'Statistics of respondents in Poland'
alignment: LRRRRRR
include: build/results/analysis/respondent-stats-POL.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.4
    - 0.2
---
```

<div class="pagebreak">

```table
---
caption: 'Statistics of respondents in Portugal'
alignment: LRRRRRR
include: build/results/analysis/respondent-stats-PRT.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.4
    - 0.2
---
```

<div class="pagebreak">

```table
---
caption: 'Statistics of respondents in Denmark'
alignment: LRRRRRR
include: build/results/analysis/respondent-stats-DNK.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.4
    - 0.2
---
```

## Model results

<div class="pagebreak">

![**Experimental design.** Probability of each pair of attribute levels to appear within the same profile. Probability is larger within attributes with fewer levels.](build/results/analysis/experimental-design.png){#fig:experimental-design}

<div class="pagebreak">

![**Posterior distributions of the varying effect of displaying options on the left-hand side compared with displaying them on the right-hand side.** (A) Expected value and uncertainty of the effect for each respondent. The dark line shows the expected value (the mean of the posterior distribution). The two shaded areas show the uncertainty (60% and 94% highest density intervals). Respondents are sorted by their expected value. (B) Posterior distribution of the population-level average.](build/results/models/hierarchical-nocovariates-nocovariances/posterior/varying/left-intercept.png){#fig:left-intercept}

<div class="pagebreak">

![**Posterior distributions of the varying effect of wind turbines.** (A) Expected value and uncertainty of the effect for each respondent. The dark line shows the expected value (the mean of the posterior distribution). The two shaded areas show the uncertainty (60% and 94% highest density intervals). Respondents are sorted by their expected value. (B) Posterior distribution of the population-level average.](build/results/models/hierarchical-nocovariates-nocovariances/posterior/varying/TECHNOLOGY___Wind.png){#fig:varying-effect-wind}

<div class="pagebreak">

![**Posterior distributions of the largest varying covariate effect across subgroups.** For each level, we show the largest covariate effect across 294 different subgroups formed by age, gender, education, and area covariates. Points show the expected value and intervals show the uncertainty (94% highest density interval).](build/results/models/hierarchical-covariates-nocovariances/subgroups/max-subgroup-effect.png){#fig:max-subgroup-effect}

<div class="pagebreak">

```table
---
caption: 'Sample statistics for parameters in the base model. {#tbl:sample-statistics}'
alignment: LRRRRR
include: build/results/models/hierarchical-nocovariates-nocovariances/posterior/diagnostics/summary-brief.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.5
    - 0.1
    - 0.1
    - 0.1
    - 0.1
    - 0.1
---
```

<div class="pagebreak">

# References
