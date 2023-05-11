# Method

## Experimental design

To assess citizen preferences for renewable electricity supply in their region, we conducted a choice experiment. Choice experiments are a widely used method to derive stated preferences in social and political science [@Hainmueller:2014; @Bansak:2021]. In our experiment, respondents are repeatedly confronted with a choice between two hypothetical designs of the electricity supply in their region. Each option is formed out of six attributes describing the electricity supply: dominant technology, land requirements, level of electricity imports into the region, household electricity prices, overhead transmission capacity expansion, and ownership of the assets (see Table XX for all attributes and attribute levels; see Figure YY for an example of the choice display). Each participant was confronted with a total of eight such choices from which we identify the relative importance of the attributes and their 25 levels. To limit biases, we randomised the combinations of attribute levels and the order in which attributes are presented (across respondents, but not within respondent).

The choice experiment was set up through conjointly.com and distributed with support from polling agencies in each country. Each respondent was exposed to eight consecutive pairs of hypothetical regional power system designs (in the following profile(s)) for a future fully renewable electricity supply with the task to select between two choices. The display of attribute levels was enhanced with small pictograms for enhancing the understanding (see Figure XX).

## Sampling

We conducted our choice experiment as an online survey between 24 January and 8 February 2022 using the platform conjointly.com. The sample was drawn from opt-in panels in the four European countries using the panel service of respondi and partners in Germany, Poland, Denmark and Portugal. We used a hard quota on the gender (50%), and a soft quota on respondents from rural areas (>=30%).

## Country selection

We chose the countries to reflect the high diversity of geographical, demographic, and socio-economic as well as historic-cultural differences across Europe55–58. The country selection is not representative for all EU countries but is meant to be illustrative of preferences in different places and contexts of the European Union.

Denmark as a small (40,000 km2) northern country with a small and highly urbanised population (5.8 million, of which 12% rural population), strong economic position (GDP: US398 billion, US68.007 GDP per capita), and high standard of living, with a far progressed energy transition (65% renewable electricity). Germany is a large (349,380 km2), central European country, with a large and urbanised population (83 million, 22% rural) and the largest economy in the EU (GDP: US4.26 trillion, GDP per capita: US51.203) and high living standard. Germany has advanced well in its energy transition (44.49% renewable electricity) and has a well-developed national energy and climate policy debate.

Portugal (GDP: US253 billion, GDP per capita: US24.567) and Poland (GDP: US679 billion, GDP per capita: US17.999) have a lower level of economic development and living standards in very different geographical locations and cultural-historic settings. Portugal is small southern EU country, with a population of 10.3 million and much advanced in the energy transition (58% renewable electricity). Poland, finally, is an eastern European country, large both in size (306,170 km2) and population (37.95 million), but with a high share of rural population (ca. 40%). The country has embarked on the energy transition pathway only recently (16% renewable power).

## Data analysis

We derive preferences from recorded choices using multinomial logit hierarchical bayes, a method commonly used for choice experiments59. It is based on the random utility theory and assumes that each option has a distinct utility to each respondent and that respondents will choose the option with higher utility to them. Within this framework, the deterministic part of the utility of an option (Vc,r) is the sum of the partworth utilities of its attribute levels (Equation 1). For example, for an option that includes wind turbines as the dominant technology in a region, the utility is the sum of the partworth utility of wind turbines and the partworth utilities of the five other attribute levels. In the analysis of choice experiments, partworth utilities are proxies for preferences: the higher the partworth utility of an attribute level, the more it is preferred.


Utility, however, is a latent variable: we did not measure it – we measured choices. To map from choices to utilities, multinomial logit uses a stochastic approach: The probability of choosing the first out of the two options depends on the deterministic part of the utilities (V1, V2 , Equation 2). In our two-option experiment, choice is modelled using a binomial distribution with the given probability.


To model heterogeneity among countries and participants, we use a nested hierarchical approach. The nested hierarchical approach allows partworth utilities to vary across countries and participants. We model heterogeneity by splitting partworth utilities in three parts: a global mean across all countries (α), a country-level deviation from the global mean (γc) using a zero-sum constraint, and a respondent-level deviation from the respective country mean (εr). The sum of these three parts forms the partworth utility for a single respondent (Equation 3). We do not model covariances between attribute levels as the additional computational complexity is restrictive.


As a diagnostic check, we capture whether respondents are more likely to choose one of the two options irrespective of the attribute levels. For example, respondents may be more likely to choose the left option in each displayed pair. We do this by adding an intercept term (I) to the first of the two options (Equation 4). As with the partworth utilities, we model this term in a hierarchical way with a population-level global mean and respondent-level deviations. We find that the intercept is very small or zero (@fig:left-intercept) suggesting that there is no impactful bias towards left or right options.


We implement our probabilistic model using PyMC [@ThomasWiecki:2023]. Code and data to reproduce our analysis are publicly available (after publication of the article). We sample the posterior distribution of all model parameters using a NUTS sampler. We run four independent chains of the Markov chain Monte Carlo to check for convergence. Each chain iterates a total of 13,000 times of which 10,000 iterations are tuning steps which we discard. Calculations were carried out on the ETH Euler cluster.


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

![**Posterior distributions of the varying effect of displaying options on the left-hand side compared with displaying them on the right-hand side.** (A) Expected value and uncertainty of the effect for each respondent. The dark line shows the expected value (the mean of the posterior distribution). The two shaded areas show the uncertainty (60% and 94% highest density intervals). Respondents are sorted by their expected value. (B) Posterior distribution of the population-level average.](build/results/models/hierarchical-nocovariates-nocovariances/posterior/varying/left-intercept.png){#fig:left-intercept}

<div class="pagebreak">

# References
