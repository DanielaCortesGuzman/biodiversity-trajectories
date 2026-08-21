# Data Documentation

## 01_Trajectories_Biodiversity slopes.csv
# Description: Contains site-level slope estimates for biodiversity metrics across 4,062 sampling sites in 23 European countries over the period 1971 to 2024

# Columns:
-Unique.ID:	unique identifier for each sampling site
-country:	country in which the sampling site is located
-latitude:	latitude of the sampling site (decimal degrees)
-longitude:	longitude of the sampling site (decimal degrees)
-Basin.ID: unique identifier of the basin containing the sampling site
-start.year:	first year of the time series
-end.year:	last year of the time series
-num.year:	number of years with observations (sampling years)
-dur.year:	duration of the time series (end.year - start.year)
-eqr.start:	predicted EQR value at the first year of the time series
-eqr.end:	predicted EQR value at the last year of the time series
-eqr.change:	change in predicted EQR over the time series (eqr.end - eqr.start)
-richness:	slope of taxonomic richness (total number of taxa)
-rich.rarefied:	slope of rarefied richness
-rich.coverage:	slope of coverage-based richness
-ab:	slope of total abundance
-fd:	slope of total functional diversity
-trait.ses:	slope of trait standardized effect sizes (SES)
-taxa.comp:	slope of taxonomic composition (turnover)
-trait.comp:	slope of trait composition (turnover)
-id.level:	taxonomic resolution used for identification
-richness_circularity: boolean; TRUE when the site's EQR calculation includes metrics potentially circular with richness
-ab_circularity: boolean; TRUE when the site's EQR calculation includes metrics potentially circular with abundance
-taxa_circularity: boolean; TRUE when the site's EQR calculation includes metrics potentially circular with taxonomic composition
-trait_circularity: boolean; TRUE when the site's EQR calculation includes metrics potentially circular with trait composition

## 02_Trajectories_gamm predictions main.csv
# Description: Contains model-predicted values of biodiversity slope responses as a function of initial degree of impact, and its total temporal degradation or recovery derived from the main analysis

# Columns:
-response: response variable
-eqr.sqrt:	squared root-transformed change in EQR
-eqr.start:	baseline EQR start values used as input for model predictions
-fit:	model-fitted values
-se.fit:	standard error of the fitted values
-ci.min:	lower bound of the 95% confidence interval for fitted values
-ci.max:	upper bound of the 95% confidence interval for fitted values

## 03_Trajectories_gamm predictions additional.csv
# Description: Contains model-predicted values of biodiversity slope responses as a function of initial degree of impact, and its total temporal degradation or recovery derived from the supplementary analysis including additional metrics

# Columns:
-response: response variable
-eqr.sqrt:	squared root-transformed change in EQR
-eqr.start:	baseline EQR start values used as input for model predictions
-fit:	model-fitted values
-se.fit:	standard error of the fitted values
-ci.min:	lower bound of the 95% confidence interval for fitted values
-ci.max:	upper bound of the 95% confidence interval for fitted values

## 04_Trajectories_gamm predictions circularity.csv
# Description: Contains model-predicted values of biodiversity slope responses as a function of initial degree of impact, and its total temporal degradation or recovery derived from the supplementary analysis assessing sensitivity to potential circularity between responses and metrics included in the EQR calculation

# Columns:
-response: response variable
-eqr.sqrt:	squared root-transformed change in EQR
-eqr.start:	baseline EQR start values used as input for model predictions
-fcircularity: boolean; TRUE when the site's EQR calculation includes metrics potentially circular with the response variable
-fit: model-fitted values
-se.fit: standard error of the fitted values
-ci.min: lower bound of the 95% confidence interval for fitted values
-ci.max: upper bound of the 95% confidence interval for fitted values

## 05_Trajectories_gamm predictions factor.csv
# Description: Contains model-predicted values of biodiversity slope responses as a function of initial degree of impact, and its total temporal degradation or recovery derived from the supplementary analysis assessing sensitivity to taxonomic resolution in taxon identification

# Columns:
-response: response variable
-eqr.sqrt:	squared root-transformed change in EQR
-eqr.start:	baseline EQR start values used as input for model predictions
-fid.level:	taxonomic resolution used for identification as a factor variable
-fit: model-fitted values
-se.fit: standard error of the fitted values
-ci.min: lower bound of the 95% confidence interval for fitted values
-ci.max: upper bound of the 95% confidence interval for fitted values

## 06_Trajectories_gamm predictions boot country.csv
# Description: Contains model-predicted values of biodiversity slope responses as a function of initial degree of impact, and its total temporal degradation or recovery derived from the supplementary analysis assessing sensitivity to the number of sampling sites per country

# Columns:
-response: response variable
-eqr.sqrt:	squared root-transformed change in EQR
-eqr.start:	baseline EQR start values used as input for model predictions
-fit_mean: mean of model-fitted values across bootstrap iterations
-se_mean: mean of standard error of the fitted values across bootstrap iterations
-ci.min: lower bound of the 95% confidence interval for fitted values
-ci.max: upper bound of the 95% confidence interval for fitted values

## 07_Trajectories_gamm predictions boot abundance.csv
# Description: Contains model-predicted values of abundance slope responses as a function of initial degree of impact, and its total temporal degradation or recovery derived from the supplementary analysis assessing sensitivity to the number of samples for degrading and recovering sites

# Columns:
-eqr.sqrt:	squared root-transformed change in EQR
-eqr.start:	baseline EQR start values used as input for model predictions
-fit_mean:	mean of model-fitted values of abundance across bootstrap iterations
-se_mean:	mean of standard error of the fitted abundance values across bootstrap iterations
-ci.min:	lower bound of the 95% confidence interval for fitted abundance calculated from the mean fitted values (fit_mean) and mean standard error (se_mean)
-ci.max:	upper bound of the 95% confidence interval for fitted abundance calculated from the mean fitted values (fit_mean) and mean standard error (se_mean)

## Trajectories_Main analysis.R
# Description: This script performs the analysis of community trajectories, including spatial autocorrelation structure and the relationship between EQR changes and trends in biodiversity metrics

## Trajectories_Analysis functions.R
# Description: This script contains the functions to perform the analysis of community trajectories in Trajectories_Main analysis.R and Trajectories_Supplementary analysis.R

## Trajectories_Supplementary analysis.R
# Description: This script performs the supplemenary analysis of community trajectories, including analyses of additional metrics, richness inclusion in EQR calculations, taxonomic resolition identification, and number of sampling site per country

## Trajectories_Plots.R
# Description: This script produces visualizations of community trajectories predicted values in relationship to EQR changes, based on the GAMM analyses performed in Trajectories_Main analysis.R and supplementary plots from Trajectories_Supplementary analysis.R

## Trajectories_Plotting functions.R
# Description: This script contains the functions to produce visualizations of community trajectories in Trajectories_Plots.R