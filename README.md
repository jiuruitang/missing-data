# missing-data
This repository contains R code for our paper - Using auxiliary marginal distributions in imputations 
for nonresponse while accounting for survey weights,
with application to estimating voter turnout.

4.2_CPS_analysis.R is the R code for the first part of Section 4, where we estimate voter turnout with 2018 CPS Data.
4.4_CPS_Measurement_Error.R is the R code for the second part of Section 4, where we incorporate measurement error.

The other R files are code for simulation study in the supplementary materials.
generate_population_logit.R explains how we generated populations for 8 simulation settings.
Chap4_designW_server.R illustrates how our framework generates imputed data sets when design weights are known.
Chap4_mix_server.R illusturates how our framework generates imputed data sets when design weights are not knonw.
Chap4_result.R explains how we get estimates and uncertainty quantifications from multiple imputed data sets.
