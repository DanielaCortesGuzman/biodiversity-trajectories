##### ----- Community trajectories analysis ----- #####

# This script performs the analysis of community trajectories,
# including spatial autocorrelation structure
# and the relationship between EQR changes and trends in biodiversity metrics

library(terra)
library(geoR)
library(spdep)
library(adespatial)
library(mgcv)
library(dplyr)
library(tidyr)


### ---Load the datasets
comm_slopes <- read.csv(file.choose(), header = TRUE, encoding = "latin1") # 01_Trajectories_Biodiversity slopes.csv

# Define factors
comm_slopes$fcountry <- factor(comm_slopes$country)
comm_slopes$fprov <- factor(comm_slopes$first.prov)
comm_slopes$fbasin <- factor(comm_slopes$Basin.ID)

# Transform EQR changes
plot(richness ~ eqr.change, data = comm_slopes, pch = 19)
pos <- which(comm_slopes$eqr.change > 0)
neg <- which(comm_slopes$eqr.change < 0)
comm_slopes$eqr.sqrt <- sign(comm_slopes$eqr.change) * sqrt(abs(comm_slopes$eqr.change))

plot(richness ~ eqr.sqrt, data = comm_slopes, pch = 19) # helps with data spread

range(comm_slopes$eqr.sqrt[pos])
range(comm_slopes$eqr.sqrt[neg])


### --- Spatial AC structure --- ###

# Project WGS84 decimal degrees onto a flat European surface (LCC Europe)
xy_wgs <- comm_slopes[, c("longitude", "latitude")]
xy_lcc <- project(as.matrix(xy_wgs), from = "epsg:4326", to = "epsg:3034")
comm_slopes$lon.lcc <- xy_lcc[, 1]
comm_slopes$lat.lcc <- xy_lcc[, 2]

# Detrended by longitude, latitude, country, and basin (already in main model)
richness_det <- as.numeric(
  resid(lm(richness ~ xy_lcc[, 1] + xy_lcc[, 2] + fcountry + fbasin, data = comm_slopes))
)

ab_det <- as.numeric(
  resid(lm(ab ~ xy_lcc[, 1] + xy_lcc[, 2] + fcountry + fbasin, data = comm_slopes))
)

taxa_det <- as.numeric(
  resid(lm(taxa.comp ~ xy_lcc[, 1] + xy_lcc[, 2] + fcountry + fbasin, data = comm_slopes))
)

trait_det <- as.numeric(
  resid(lm(trait.comp ~ xy_lcc[, 1] + xy_lcc[, 2] + fcountry + fbasin, data = comm_slopes))
)


# Variograms
vgram <- variog(coords = xy_lcc, data = scale(richness_det))
plot(vgram$v ~ vgram$u, type = "b", pch = 19, cex = 2)

vgram <- variog(coords = xy_lcc, data = scale(ab_det))
plot(vgram$v ~ vgram$u, type = "b", pch = 19, cex = 2)

vgram <- variog(coords = xy_lcc, data = scale(taxa_det))
plot(vgram$v ~ vgram$u, type = "b", pch = 19, cex = 2)

vgram <- variog(coords = xy_lcc, data = scale(bio_det))
plot(vgram$v ~ vgram$u, type = "b", pch = 19, cex = 2)

# MEMs based on fine-scale thresholds
thresh1 <- 0
thresh2 <- c(100, 500, 1000, 5000, 10000, 25000, 50000, 100000)

can_dnear <- lapply(thresh2, function(d2) {
  nb <- dnearneigh(xy_lcc, d1 = 0, d2 = d2)
  nb2listw(nb, style = "B", zero.policy = TRUE)
})
names(can_dnear) <- thresh2

# MEM selection: >1% R2 and significant
(dnear_richness <- listw.select(
  richness_det, can_dnear,
  MEM.autocor = "positive", p.adjust = TRUE
))
richness_best <- as.data.frame(
  dnear_richness$best$MEM.select[, c(1:94)]
) # no MEM explained >1%

(dnear_ab <- listw.select(
  ab_det, can_dnear,
  MEM.autocor = "positive", p.adjust = TRUE
))
ab_best <- as.data.frame(
  dnear_ab$best$MEM.select[, c(1:92)]
) # no mem explained >1%

(dnear_taxa <- listw.select(
  taxa_det, can_dnear,
  MEM.autocor = "positive", p.adjust = TRUE
))
taxa_best <- as.data.frame(
  dnear_taxa$best$MEM.select[, c(1:46)]
) # no mem explained >1%

(dnear_trait <- listw.select(
  trait_det, can_dnear,
  MEM.autocor = "positive", p.adjust = TRUE
))
trait_best <- as.data.frame(
  dnear_trait$best$MEM.select[, c(1:40)]
) # no mem explained >1%


### --- GAMMs --- ###

# Load fitting and predicting functions
source("R Scripts/Trajectories_Analysis functions.R")

## -- Richness
# Fit
richness_gamm <- fit_gamm(
  response = "richness",
  data = comm_slopes,
  dist_family = "gaussian"
)
gam.check(richness_gamm$gam)
summary(richness_gamm$gam)

# Predict
richness_rec <- predict_gamm(
  model = richness_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos
)

richness_deg <- predict_gamm(
  model = richness_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg
)
richness_predicted <- merge(richness_rec, richness_deg, all = TRUE)



## -- Abundance
# Fit
ab_gamm <- fit_gamm(
  response = "ab",
  data = comm_slopes,
  dist_family = "gaussian"
)
gam.check(ab_gamm$gam)
summary(ab_gamm$gam)

# Predict
ab_rec <- predict_gamm(
  model = ab_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos
)

ab_deg <- predict_gamm(
  model = ab_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg
)
ab_predicted <- merge(ab_rec, ab_deg, all = TRUE)



# Functional diversity
# Fit
fd_gamm <- fit_gamm(
  response = "fd",
  data = comm_slopes,
  dist_family = "gaussian"
)
gam.check(fd_gamm$gam)
summary(fd_gamm$gam)

# Predict
fd_rec <- predict_gamm(
  model = fd_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos
)

fd_deg <- predict_gamm(
  model = fd_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg
)
fd_predicted <- merge(fd_rec, fd_deg, all = TRUE)



# Trait SES
# Fit
trait_ses_gamm <- fit_gamm(
  response = "trait.ses",
  data = comm_slopes,
  dist_family = "gaussian"
)
gam.check(trait_ses_gamm$gam)
summary(trait_ses_gamm$gam)

# Predict
ses_rec <- predict_gamm(
  model = trait_ses_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos
)

ses_deg <- predict_gamm(
  model = trait_ses_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg
)
ses_predicted <- merge(ses_rec, ses_deg, all = TRUE)


# Taxonomic composition
# Fit
taxa_comp_gamm <- fit_gamm(
  response = "taxa.comp",
  data = comm_slopes,
  dist_family = "gaussian"
)
gam.check(taxa_comp_gamm$gam)
summary(taxa_comp_gamm$gam)

# Predict
taxa_comp_rec <- predict_gamm(
  model = taxa_comp_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos
)

taxa_comp_deg <- predict_gamm(
  model = taxa_comp_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg
)
taxa_comp_predicted <- merge(taxa_comp_rec, taxa_comp_deg, all = TRUE)



# Trait composition
# Fit
trait_comp_gamm <- fit_gamm(
  response = "trait.comp",
  data = comm_slopes,
  dist_family = gaussian
)
gam.check(trait_comp_gamm$gam)
summary(trait_comp_gamm$gam)

# Predict
trait_comp_rec <- predict_gamm(
  model = trait_comp_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos
)

trait_comp_deg <- predict_gamm(
  model = trait_comp_gamm,
  data = comm_slopes,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg
)
trait_comp_predicted <- merge(trait_comp_rec, trait_comp_deg, all = TRUE)



### --- Export predictions --- ###
# Combine all predictions into a single dataframe for export
dfs <- list(
  richness    = richness_predicted,
  ab    = ab_predicted,
  fd   = fd_predicted,
  ses   = ses_predicted,
  taxa  = taxa_comp_predicted,
  trait = trait_comp_predicted
)
dfs <- purrr::map2(
  dfs,
  names(dfs),
  ~ mutate(.x, response = .y)
)
all_predictions <- bind_rows(dfs)
all_predictions <- all_predictions[, -c(3:8)]
all_predictions <- all_predictions |> relocate(response)

# Export
write.csv(all_predictions, "02_Trajectories_gamm predictions main.csv", row.names = FALSE)


### --- Checking points: influence of fixed covariates --- ###
## -- Richness
# Fit
richness_gamm_fixed <- fit_gamm_fixed(
  response = "richness",
  data = comm_slopes,
  dist_family = "gaussian"
)
# Partial R2
summary(richness_gamm$gam)$r.sq - summary(richness_gamm_fixed$gam)$r.sq


## -- Abundance
# Fit
ab_gamm_fixed <- fit_gamm_fixed(
  response = "ab",
  data = comm_slopes,
  dist_family = "gaussian"
)
# Partial R2
summary(ab_gamm$gam)$r.sq - summary(ab_gamm_fixed$gam)$r.sq


# Functional diversity
# Fit
fd_gamm_fixed <- fit_gamm_fixed(
  response = "fd",
  data = comm_slopes,
  dist_family = "gaussian"
)
# Partial R2
summary(fd_gamm$gam)$r.sq - summary(fd_gamm_fixed$gam)$r.sq


# Trait SES
# Fit
trait_ses_gamm_fixed <- fit_gamm_fixed(
  response = "trait.ses",
  data = comm_slopes,
  dist_family = "gaussian"
)
# Partial R2
summary(trait_ses_gamm$gam)$r.sq - summary(trait_ses_gamm_fixed$gam)$r.sq


# Taxonomic composition
# Fit
taxa_comp_gamm_fixed <- fit_gamm_fixed(
  response = "taxa.comp",
  data = comm_slopes,
  dist_family = "gaussian"
)
# Partial R2
summary(taxa_comp_gamm$gam)$r.sq - summary(taxa_comp_gamm_fixed$gam)$r.sq


# Trait composition
# Fit
trait_comp_gamm_fixed <- fit_gamm_fixed(
  response = "trait.comp",
  data = comm_slopes,
  dist_family = "gaussian"
)
# Partial R2
summary(trait_comp_gamm$gam)$r.sq - summary(trait_comp_gamm_fixed$gam)$r.sq
