# Introduction

The following is a brief analysis of our choice experiment results with ~1000 respondents from each of the four countries: Germany, Poland, Denmark, Portugal.

# Hypotheses

* 🚫 H1: People prefer high levels of local self-sufficiency even if they must pay a premium. @sec:H1
* 🚫 H2: People value ownership higher than prices. @sec:H2
* 🚫 H3: People value land requirements higher than prices. @sec:H3
* ✅/🚫 H4: People prefer low levels of land requirements even if they must pay a premium. @sec:H4
* ✅ H5: People are willing to pay a premium for shared ownership. @sec:H5
* 🚫 H6: People in urban areas show a higher preference for land requirements than people from rural areas. @sec:H6
* ✅/🚫 H7: People prefer high levels of self-sufficiency even if they must pay a land premium. @sec:H7
* 🚫 H8: People prefer land requirements more if dominant technology is solar. @sec:H8
* ✅ H9: People prefer imports more if dominant technology is wind. @sec:H9
* H10: @sec:H10
* ✅ H11: People in wealthy countries show a stronger preference for local generation infrastructure. @sec:H11
* ✅ H12: People in countries with high deployment levels of renewables show a stronger preference for local generation infrastructure. @sec:H12
* 🚫 H13: People in countries with low population density show a stronger preference for local generation infrastructure. @sec:H13
* 🚫 H14: People in resource-rich countries show a stronger preference for local generation infrastructure. @sec:H14
* H15: @sec:H15

## H1: People prefer high levels of local self-sufficiency even if they must pay a premium. {#sec:H1}

No.

I create a statistical model that contains an interaction between PRICES and SHARE_IMPORTS [@Hainmueller:2014]. If this hypothesis was correct, the Average Component Interaction Effect (ACIE) of high prices would need to increase with higher self-sufficiency rates. This is not the case (bottom twelve lines in @fig:H1).

![**AMCEs and ACIE for an interaction between PRICES and SHARE_IMPORTS.**](build/H1.png){#fig:H1}

## H2: People value ownership higher than prices. {#sec:H2}

There is no formal way to test this hypothesis. Likely, it does not hold.

We could test this hypothesis by comparing the magnitude of the effect sizes of the attributes. Doing so, ownership's largest effect size is ~5 percentage points (from public to private), and prices' largest effect size is ~3 percentage points (from +0% to +60%) (@fig:H2). We must not forget that these effect sizes are somewhat arbitrary, as they are based on our arbitrary choice of levels. Still, as the effecct sizes differ quite a lot, we can likely reject this hypothesis.

![**AMCEs for all levels.**](build/H2.png){#fig:H2}

## H3: People value land requirements higher than prices. {#sec:H3}

There is no formal way to test this hypothesis. Likely, it does not hold.

We could test this hypothesis by comparing the magnitude of the effect sizes of the attributes. Doing so, land requirements's largest effect size is <5 percentage points, and prices' largest effect size is ~30 percentage points (from +0% to +60%) (@fig:H2). We must not forget that these effect sizes are somewhat arbitrary, as they are based on our arbitrary choice of levels. Still, as the effecct sizes differ quite a lot, we can likely reject this hypothesis.

## H4: People prefer low levels of land requirements even if they must pay a premium. {#sec:H4}

Difficult case. The effect likely exists, but is difficult to prove.

I create a statistical model that contains an interaction between PRICES and LAND [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions involving high prices (+30%, +45%, +60%) reduces with higher land requirements (bottom 12 lines in @fig:H4) suggesting that people are willing to accept higher prices if land requirements are low. The effect is borderline statistically significant.

The interaction effect between PRICES and LAND are generally extremely high -- up to 30 percentage points. I do not understand at this point why. We need to find out.

![**AMCEs and ACIE for an interaction between PRICES and LAND.**](build/H4.png){#fig:H4}

## H5: People are willing to pay a premium for shared ownership. {#sec:H5}

Yes. People are willing to pay a premium for public and community compared to private ownership. This effect is strong and statistically significant.

I create a statistical model that contains an interaction between PRICES and OWNERSHIP [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions involving private ownership are negative (four of the eight bottom lines in @fig:H5 -- sorry I cannot improve the confusing visualisation at this point). This suggests that people are willing to pay a cost premium for shared (public / community) ownership. With ~10 percentage points higher probability for shared ownership, this effect is strong and statistically significant.

![**AMCEs and ACIE for an interaction between PRICES and OWNERSHIP.**](build/H5.png){#fig:H5}

## H6: People in urban areas show a higher preference for land requirements than people from rural areas. {#sec:H6}

The data do not support this hypothesis.

I am using the method to measure subgroup preferences described by @Leeper:2020. I remove all respondents that gave no answer to the question about their area. I then calculate differences in marginal means between the urban and the rural population for all attribute levels. All differences including the one about land requirements are very small (green lines in @fig:H6) and only for two levels (unrelated to land requirements) can we reject the null hypothesis that there is an effect.

![**Differences in marginal means between urban and rural population.**](build/H6.png){#fig:H6}

## H7: People prefer high levels of self-sufficiency even if they must pay a land premium. {#sec:H7}

Borderline.

I create a statistical model that contains an interaction between SHARE_IMPORTS and LAND [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of high self sufficiency (import share 10%) are generally high (line 9--12 from the bottom of @fig:H7), suggesting that people favor self-sufficiency. The ACIE of very high land requirements is about 5% less probable than of all other land requirements suffesting that people disfavor very high land requirements even for high self-sufficiency. The effect is not strong though and not statistically significant.

![**AMCEs and ACIE for an interaction between SHARE_IMPORTS and LAND.**](build/H7.png){#fig:H7}

## H8: People prefer land requirements more if dominant technology is solar. {#sec:H8}

No. The opposite, if at all.

I create a statistical model that contains an interaction between LAND and TECHNOLOGY [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions containing open-field PV are the lowest (lines 5--8 from the bottom of @fig:H8) suggesting that people disfavor land requirements even more when open-field PV is dominant. Wind has a similar interaction with land requirements, albeit generally weaker. The interaction effect is generally strong, with up to 10% lower probability of people choosing higher land requirements because of the dominant technnology.

This result is confusing. It may be an artefact: AMCE's of land requirements are generally very small, so it may as well be that this can be ignored. However, we need to find an explanation (similar phenomena as in @sec:H4).

It is possible that respondents associated land competition with agriculture with open-field PV, leading to the surprisingly poor results of open-field PV. We did not exclude land competition.

![**AMCEs and ACIE for an interaction between LAND and TECHNOLOGY.**](build/H8.png){#fig:H8}

## H9: People prefer imports more if dominant technology is wind. {#sec:H9}

Yes.

I create a statistical model that contains an interaction between SHARE_IMPORTS and TECHNOLOGY [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions containing wind are positive (bottom three lines in @fig:H9) suggesting that people choose higher import shares when wind is in the profile. Wind increases the probability of higher import shares being chosen by ~5--10 percentage points, and this effect is statistically significant on the 5% level.

![**AMCEs and ACIE for an interaction between SHARE_IMPORTS and TECHNOLOGY.**](build/H9.png){#fig:H9}

## H10 {#sec:H10}

## H11: People in wealthy countries show a stronger preference for local generation infrastructure. {#sec:H11}

The data support this hypothesis.

I am using the method to measure subgroup preferences described by @Leeper:2020. I then calculate conditional marginal means for all attribute levels.

Country GDP per capita leads to the following order of countries (Worldbank data):

1. Denmark ($56000)
2. Germany ($41000)
3. Portugal ($20000)
4. Poland ($15000)

Germany has the highest preference for no imports, followed by Denmark, Portugal, and Poland (@fig:H11). Some of these differences are statistically significant.

Germany has the lowest preference for 90% imports, followed by Denmark, Portugal, and Poland (@fig:H11). Some of these differences are statistically significant.

![**Conditional marginal means for all countries.**](build/H11.png){#fig:H11}

## H12: People in countries with high deployment levels of renewables show a stronger preference for local generation infrastructure. {#sec:H12}

The data support this hypothesis.

We define deployment levels as the ratio of 2019 deployment of solar and wind renewables [@IRENA:2020a] to electricity demand [@Trondle:2019].

Deployment:

* Denmark  6117 +    1079 =   7196 MW
* Germany  60840 +  49018 = 109858 MW
* Poland   5917 +    1300 =   7217 MW
* Portugal 5233 +     842 =   6075 MW

Demand:

* Germany  (493 TWh yr^-1^)
* Denmark  ( 33 TWh yr^-1^)
* Poland   (168 TWh yr^-1^)
* Portugal ( 50 TWh yr^-1^)

Leads to the following country ranking:

1. Germany  223 MW TWh^-1^ yr^1^
2. Denmark  218 MW TWh^-1^ yr^1^
3. Portugal 122 MW TWh^-1^ yr^1^
4. Poland    43 MW TWh^-1^ yr^1^

Germany has highest deployment and highest preferences, followed by Denmark, Portugal, and Poland.
This ranking is almos the same as in H11, and therefore the data support the hypothesis.

## H13: People in countries with low population density show a stronger preference for local generation infrastructure. {#sec:H13}

The data do not support this hypothesis.

Population density of the four countries (Wikipedia):

1. Portugal (112 km^-2^)
2. Poland (123 km^-2^)
3. Denmark (138 km^-2^)
4. Germany (232 km^-2^)

This is in analogy to H11, however the country ranking is almost completely inverted. The conditional marginal means indicate the exact opposite (@fig:H11).

## H14: People in resource-rich countries show a stronger preference for local generation infrastructure. {#sec:H14}

The data do not support this hypothesis.

We define resource-richness as the ratio of generation potential to electricity demand. Data taken from @Trondle:2019. The method is in analogy to H11.

* Germany  (1 183 and 493 TWh yr^-1^)
* Denmark  (  509 and  33 TWh yr^-1^)
* Poland   (  548 and 168 TWh yr^-1^)
* Portugal (  363 and  50 TWh yr^-1^)

Leads to the following ranking of resource-richness:

1. Denmark  (15.4)
2. Portugal ( 7.2)
3. Poland   ( 3.3)
4. Germany  ( 2.4)

Germany shows higher preference for local generation infrastructure than Portugal and Poland, but has fewer resources (@fig:H11). Denmark has by far the best resources and generally rather high preferences, but the two are not in proportion.

## H15 {#sec:H15}

# Appendix

Additional tables and figures in arbitrary order.

## Respondents

```table
---
caption: 'Statistics of respondents in Germany'
alignment: LRRRRRR
include: build/DEU/respondent-stats.csv
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
include: build/POL/respondent-stats.csv
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
include: build/PRT/respondent-stats.csv
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
include: build/DNK/respondent-stats.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.4
    - 0.2
---
```

<div class="pagebreak">

## AMCE

![**AMCE values of all levels for Germany.**](build/DEU/amce.png){#fig:amce-deu}

![**AMCE values of all levels for Poland.**](build/POL/amce.png){#fig:amce-pol}

![**AMCE values of all levels for Portugal.**](build/PRT/amce.png){#fig:amce-por}

![**AMCE values of all levels for Denmark.**](build/DNK/amce.png){#fig:amce-den}

![**AMCE values of all levels for all four countries.**](build/amce.png){#fig:amce}

## Marginal Means

![**Marginal Mean values of all levels for Germany.**](build/DEU/mm.png){#fig:mm-deu}

![**Marginal Mean values of all levels for Poland.**](build/POL/mm.png){#fig:mm-pol}

![**Marginal Mean values of all levels for Portugal.**](build/PRT/mm.png){#fig:mm-por}

![**Marginal Mean values of all levels for Denmark.**](build/DNK/mm.png){#fig:mm-den}

![**Marginal Mean values of all levels for all four countries.**](build/mm.png){#fig:mm}

# References
