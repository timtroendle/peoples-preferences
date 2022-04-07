# Introduction

The following is a brief analysis of our choice experiment results with ~1000 respondents from each of the four countries: Germany, Poland, Denmark, Portugal.

# Hypotheses

The data support some hypothesis (✅) and do not support others (🚫). We cannot (right now) test some hypothesis (❓).

* 🚫 H1: People prefer high levels of local self-sufficiency even if they must pay a premium. @sec:H1
* 🚫 H2: People value ownership higher than prices. @sec:H2
* 🚫 H3: People value land requirements higher than prices. @sec:H3
* ✅/🚫 H4: People prefer low levels of land requirements even if they must pay a premium. @sec:H4
* ✅ H5: People are willing to pay a premium for shared ownership. @sec:H5
* 🚫 H6: People in urban areas show a higher preference for land requirements than people from rural areas. @sec:H6
* ✅/🚫 H7: People prefer high levels of self-sufficiency even if they must pay a land premium. @sec:H7
* 🚫 H8: People prefer land requirements more if dominant technology is solar. @sec:H8
* ✅ H9: People prefer imports more if dominant technology is wind. @sec:H9
* ❓ H10: There are no differences between countries but differences in regards to demographic and regional conditions. @sec:H10
* ✅ H11: People in wealthy countries show a stronger preference for local generation infrastructure. @sec:H11
* ✅ H12: People in countries with high deployment levels of renewables show a stronger preference for local generation infrastructure. @sec:H12
* 🚫 H13: People in countries with low population density show a stronger preference for local generation infrastructure. @sec:H13
* 🚫 H14: People in resource-rich countries show a stronger preference for local generation infrastructure. @sec:H14
* ✅/❓H15: Preference for self-sufficiency varies between countries (this effect is not explained by demographics and regional differences). @sec:H15
* 🚫 H16: Gender interacts with preferences. @sec:H16
    * 🚫 H16a: Men prefer lower land requirements than women.
    * 🚫 H16b: Men prefer lower prices than women.
    * ❓/🚫 H16c: Men have no technology preference, whereas women prefer solar.
    * 🚫 H16d: Women prefer lower imports than men.
    * 🚫 H16e: Women prefer, more clearly than men, public over private ownership.
* ✅ H17: Age interacts with preferences. @sec:H17
    * 🚫 H17a: People over 55 (old people) prefer lower land requirements than people below 30 (young people).
    * 🚫 H17b: People over 55 prefer lower transmission infrastructure than people below 30.
    * ✅ H17c: People below 30 prefer high imports.
    * ✅ H17d: People over 55 prefer low imports.
    * ✅ H17e: There is no difference in technology preferences across ages.
    * 🚫 H17f: People over 55 prefer higher transmission than younger people.
* 🚫 H18: Place attachment interacts with preferences. @sec:H18
    * 🚫 H18a: People prefer lower land requirements/ transmission infrastructure the longer they live in the region. = People with strong place attachment (+15 years in the region) prefer lower land requirements and lower transmission infrastructure than people with “weaker” place attachment (1-5 years in the region).
    * 🚫 H18b: People with strong place attachment prefer higher imports.
    * 🚫 H18c: People living more than 15 years in the region prefer solar technologies.
* ✅ H19: Income interacts with preferences. @sec:H19
    * 🚫 H19a: People with high income prefer community ownership more than people with lower income.
    * 🚫 H19b: People with low income prefer lower land requirements than people with higher income.
    * 🚫 H19c: People with low income prefer lower prices than people with higher income.
* 🚫 H20: Education interacts with preferences. @sec:H20
    * 🚫 H20a: People with lower education prefer lower imports, than people with higher education.
    * 🚫 H20b: People with higher education prefer higher land use requirements (understand the necessity to act)
* ✅ H21: Climate concern interacts with preferences. @sec:H21
    * 🚫 People with high climate concerns prefer higher land requirements and transmission infrastructure.
* H22: Partisanship interacts with preferences. @sec:H22
    * H22a: People supporting right wing parties prefer lower land requirements. (don’t support the energy transition). But, prefer lower imports.
    * H22b: People supporting green parties prefer higher land requirements.
* ✅ H23: Current deployment interacts with preferences. @sec:H23
    * 🚫 H23a: People prefer (higher) land requirements more if they have renewables in their region (wind, solar or both)/ People prefer lower land requirements if they have no renewables in their region.
    * ✅ H23b: People prefer technology wind, if they have wind in their region.
    * ✅ H23c: People prefer technology open space solar (Freiflächenanlagen), if they have open space solar in their region.
    * 🚫 H23d: People prefer (expansion) transmission infrastructure more if they have no renewables in their region.
* ✅ H24: People prefer higher land requirements more if ownership is public. @sec:H24
* ✅ H25: People prefer higher imports more if transmission is high. @sec:H25
* 🚫 H26: People prefer higher land requirements more if transmission infrastructure is low. @sec:26
* 🚫/✅ H27: People prefer higher transmission infrastructure more when prices are low. @sec:27
* 🚫/✅ H28: People prefer higher transmission infrastructure more when the dominant technology is wind. @sec:28
* 🚫 H29: Technology preferences are independent of ownership. @sec:29


## H1: People prefer high levels of local self-sufficiency even if they must pay a premium. {#sec:H1}

No.

I create a statistical model that contains an interaction between PRICES and SHARE_IMPORTS [@Hainmueller:2014]. If this hypothesis was correct, the Average Component Interaction Effect (ACIE) of high prices would need to increase with higher self-sufficiency rates. This is not the case (bottom twelve lines in @fig:H1).

![**AMCEs and ACIEs for an interaction between PRICES and SHARE_IMPORTS.**](build/H1.png){#fig:H1}

## H2: People value ownership higher than prices. {#sec:H2}

There is no formal way to test this hypothesis. Likely, it does not hold.

We could test this hypothesis by comparing the magnitude of the effect sizes of the attributes. Doing so, ownership's largest effect size is ~5 percentage points (from public to private), and prices' largest effect size is ~3 percentage points (from +0% to +60%) (@fig:H2). We must not forget that these effect sizes are somewhat arbitrary, as they are based on our arbitrary choice of levels. Still, as the effecct sizes differ quite a lot, we can likely reject this hypothesis.

![**AMCEs for all levels.**](build/H2.png){#fig:H2}

## H3: People value land requirements higher than prices. {#sec:H3}

There is no formal way to test this hypothesis. Likely, it does not hold.

We could test this hypothesis by comparing the magnitude of the effect sizes of the attributes. Doing so, land requirements's largest effect size is <5 percentage points, and prices' largest effect size is ~30 percentage points (from +0% to +60%) (@fig:H2). We must not forget that these effect sizes are somewhat arbitrary, as they are based on our arbitrary choice of levels. Still, as the effecct sizes differ quite a lot, we can likely reject this hypothesis.

## H4: People prefer low levels of land requirements even if they must pay a premium. {#sec:H4}

The effect likely exists, but is too small to prove.

I create a statistical model that contains an interaction between PRICES and LAND [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions involving high prices (+30%, +45%, +60%) reduces with higher land requirements (bottom 12 lines in @fig:H4) suggesting that people are willing to accept higher prices if land requirements are low. The effect is borderline statistically significant.

The interaction effect between PRICES and LAND are generally extremely high -- up to 30 percentage points. I do not understand at this point why. We need to find out.

![**AMCEs and ACIEs for an interaction between PRICES and LAND.**](build/H4.png){#fig:H4}

## H5: People are willing to pay a premium for shared ownership. {#sec:H5}

Yes. People are willing to pay a premium for public and community compared to private ownership. This effect is strong and statistically significant.

I create a statistical model that contains an interaction between PRICES and OWNERSHIP [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions involving private ownership are negative (four of the eight bottom lines in @fig:H5 -- sorry I cannot improve the confusing visualisation at this point). This suggests that people are willing to pay a cost premium for shared (public / community) ownership. With ~10 percentage points higher probability for shared ownership, this effect is strong and statistically significant.

![**AMCEs and ACIEs for an interaction between PRICES and OWNERSHIP.**](build/H5.png){#fig:H5}

## H6: People in urban areas show a higher preference for land requirements than people from rural areas. {#sec:H6}

The data do not support this hypothesis.

I am using the method to measure subgroup preferences described by @Leeper:2020. I remove all respondents that gave no answer to the question about their area. I then calculate differences in marginal means between the urban and the rural population for all attribute levels. All differences including the one about land requirements are very small (green lines in @fig:H6) and only for two levels (unrelated to land requirements) can we reject the null hypothesis that there is an effect.

![**Differences in marginal means between rural and urban population.**](build/H6.png){#fig:H6}

## H7: People prefer high levels of self-sufficiency even if they must pay a land premium. {#sec:H7}

Borderline.

I create a statistical model that contains an interaction between SHARE_IMPORTS and LAND [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of high self sufficiency (import share 10%) are generally high (line 9--12 from the bottom of @fig:H7), suggesting that people favor self-sufficiency. The ACIE of very high land requirements is about 5% less probable than of all other land requirements suffesting that people disfavor very high land requirements even for high self-sufficiency. The effect is not strong though and not statistically significant.

![**AMCEs and ACIEs for an interaction between SHARE_IMPORTS and LAND.**](build/H7.png){#fig:H7}

## H8: People prefer land requirements more if dominant technology is solar. {#sec:H8}

No. The opposite, if at all.

I create a statistical model that contains an interaction between LAND and TECHNOLOGY [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions containing open-field PV are the lowest (lines 5--8 from the bottom of @fig:H8) suggesting that people disfavor land requirements even more when open-field PV is dominant. Wind has a similar interaction with land requirements, albeit generally weaker. The interaction effect is generally strong, with up to 10% lower probability of people choosing higher land requirements because of the dominant technnology.

This result is confusing. It may be an artefact: AMCE's of land requirements are generally very small, so it may as well be that this can be ignored. However, we need to find an explanation (similar phenomena as in @sec:H4).

It is possible that respondents associated land competition with agriculture to open-field PV, leading to the surprisingly poor results of open-field PV. We did not exclude land competition.

![**AMCEs and ACIEs for an interaction between LAND and TECHNOLOGY.**](build/H8.png){#fig:H8}

## H9: People prefer imports more if dominant technology is wind. {#sec:H9}

Yes.

I create a statistical model that contains an interaction between SHARE_IMPORTS and TECHNOLOGY [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions containing wind are positive (bottom three lines in @fig:H9) suggesting that people choose higher import shares when wind is in the profile. Wind increases the probability of higher import shares being chosen by ~5--10 percentage points, and this effect is statistically significant on the 5% level.

![**AMCEs and ACIEs for an interaction between SHARE_IMPORTS and TECHNOLOGY.**](build/H9.png){#fig:H9}

## H10: There are no differences between countries but differences in regards to demographic and regional conditions. {#sec:H10}

We will not be able to test this. We cannot show causality for non-treatment data like demographics.

However, by using other models than described in [@Hainmueller:2014; @Leeper:2020], we can can get some evidence whether this may be true. For that, we need to build two models; both including demographics but only one including the country as explaining variable. We can then test whether the model including country performs better than the one without. If it does, we have evidence that this hypothesis can be rejected.

Less formally, I find it striking how similar the preferences in the four countries are in all attributes but PRICES and SHARE_IMPORTS. Even odd effects in LAND and TRANSMISSION that we see in the pooled data exist in all countries (@fig:H9). We should understand this better.

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

## H15: Preference for self-sufficiency varies between countries (this effect is not explained by demographics and regional differences). {#sec:H15}

Yes, preferences for self-sufficiency vary between countries. The effect is strongest for the extreme levels of self-sufficienncy and the difference is largest for Germany and Poland (@fig:H2). While German respondents favor high self-sufficiency most (marginal mean of no imports > 60%) and favor low self-sufficiency least (marginal mean of 90% imports < 40%), Polish respondents are more indifferent towards self-sufficiency (marginal means of 55% and 45%, respectively).

Right now, and with the tools we have [@Hainmueller:2014; @Leeper:2020], we cannot show that this effect is not explained best by demographics and regional differences. We need a logit model for that (as typically used in market research) (see also @sec:H10).

## H16: Gender interacts with preferences. {#sec:H16}

No.

I am using the method to measure subgroup preferences described by @Leeper:2020. I remove all respondents that answered "other". I then calculate differences in marginal means between the male and the female respondents for all attribute levels. All differences are very small (@fig:H16) and almost all are not statisticall significant.

![**Differences in marginal means between male and female respondents.**](build/H16.png){#fig:H16}

## H17: Age interacts with preferences. {#sec:H17}

Yes.

I am using the method to measure subgroup preferences described by @Leeper:2020. I then calculate conditional marginal means for all attribute levels.

Old people prefer lower imports and lower prices more than younger people. They also like public ownership and PV more then younger, but these are not statistically significant.

![**Marginal means conditional to age of respondents.**](build/H17.png){#fig:H17}

## H18: Place attachment interacts with preferences. {#sec:H18}

No.

I am using the method to measure subgroup preferences described by @Leeper:2020. I then calculate conditional marginal means for all attribute levels.

![**Marginal means conditional to place attachment of respondents.**](build/H18.png){#fig:H18}

## H19: Income interacts with preferences. {#sec:H19}

Yes.

I am using the method to measure subgroup preferences described by @Leeper:2020. I then calculate conditional marginal means for all attribute levels.

People with higher income do actually prefer community ownership less, but they do prefer public ownership more then people with lower income (@fig:H19). This effect is statistically not significant.

There is no clear and certainly no statistically significant interaction with land requirements (@fig:H19).

People with low income do _not_ prefer lower prices more than people with high income. BUT: People with low income prefer higher prices more than people with high income. That means that medium and high income respondents are more price-sensitive.

The _high_ income group is very small leading to very large confidence intervals.

![**Marginal means conditional to income of respondents.**](build/H19.png){#fig:H19}

## H20: Education interacts with preferences. {#sec:H20}

No.

I am using the method to measure subgroup preferences described by @Leeper:2020. I then calculate conditional marginal means for all attribute levels.

![**Marginal means conditional to education of respondents.**](build/H20.png){#fig:H20}

## H21: Climate concern interacts with preferences. {#sec:H21}

Yes.

I am using the method to measure subgroup preferences described by @Leeper:2020. I then calculate conditional marginal means for all attribute levels.

Respondents with high climate concern do _not_ prefer land requirements of transmission more than respondents with low climate concern, but they _do_ prefer self sufficiency more (@fig:H21).

![**Marginal means conditional to climate concern of respondents.**](build/H21.png){#fig:H21}

## H22: Partisanship interacts with preferences. {#sec:H22}

## H23: Current deployment interacts with preferences. {#sec:H23}

Yes.

I am using the method to measure subgroup preferences described by @Leeper:2020. I then calculate conditional marginal means for all attribute levels.

Current deployment does not interact with land requirements or transmission infrastructure (@fig:H23). But people with wind turbines in their region prefer wind turbines more. People with open-field PV prefer open-field PV more. The effect is stronger for wind than for PV.

![**Marginal means conditional to current renewable deployment in the region of respondents.**](build/H23.png){#fig:H23}

## H24: People prefer higher land requirements more if ownership is public. {#sec:H24}

Yes, but the effect is small.

I create a statistical model that contains an interaction between LAND and OWNERSHIP [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of all interactions containing private ownership and high land requirements are negative (bottom three lines in @fig:H9) suggesting that people choose higher land requirements when ownership is public. The effect is small. Community ownership is worse than private which is worse than public.

![**AMCEs and ACIEs for an interaction between LAND and OWNERSHIP.**](build/H24.png){#fig:H24}

## H25: People prefer higher imports more if transmission is high. {#sec:H25}

Yes.

I create a statistical model that contains an interaction between SHARE_IMPORTS and TRANSMISSION [@Hainmueller:2014]. The Average Component Interaction Effect (ACIE) of interactions containing high import shares (50%, 90%) are higher for large transmission expansion (75%) (@fig:H9) suggesting that people choose higher import shares when transmission is high (and vice versa).

![**AMCEs and ACIEs for an interaction between SHARE_IMPORTS and TRANSMISSION.**](build/H25.png){#fig:H25}

## H26: People prefer higher land requirements more if transmission infrastructure is low. {#sec:26}

No clear effect.

I create a statistical model that contains an interaction between LAND and TRANSMISSION [@Hainmueller:2014].

![**AMCEs and ACIEs for an interaction between LAND and TRANSMISSION.**](build/H26.png){#fig:H26}

## H27: People prefer higher transmission infrastructure more when prices are low. {#sec:27}

Yes, but the effect is not very clear.

I create a statistical model that contains an interaction between TRANSMISSION and PRICE [@Hainmueller:2014].

![**AMCEs and ACIEs for an interaction between TRANSMISSION and PRICE.**](build/H27.png){#fig:H27}

## H28: People prefer higher transmission infrastructure more when the dominant technology is wind. {#sec:28}

Maybe a small effect.

I create a statistical model that contains an interaction between TRANSMISSION and TECHNOLOGY [@Hainmueller:2014].

![**AMCEs and ACIEs for an interaction between TRANSMISSION and TECHNOLOGY.**](build/H28.png){#fig:H28}

## H29: Technology preferences are independent of ownership. {#sec:29}

No, there is an effect.

I create a statistical model that contains an interaction between TRANSMISSION and TECHNOLOGY [@Hainmueller:2014]. Open-field PV is preferred less when ownerhsip is community or private.

![**AMCEs and ACIEs for an interaction between TECHNOLOGY and OWNERSHIP.**](build/H29.png){#fig:H29}

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
