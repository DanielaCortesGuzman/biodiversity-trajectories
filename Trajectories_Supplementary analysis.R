##### ----- Community trajectories analysis ----- #####

# This script performs the supplementary analysis of community trajectories,
# including analyses of additional metrics,
# and sensitivity analyses to test the robustness of the results

library(mgcv)
library(dplyr)
library(tidyr)
library(terra)
library(gratia)

### ---Load the datasets
comm_slopes <- read.csv(file.choose(), header = TRUE, encoding = "latin1") # 01_Trajectories_biodiversity slopes.csv

# Define factors
comm_slopes$fcountry <- factor(comm_slopes$country)
comm_slopes$fprov <- factor(comm_slopes$first.prov)
comm_slopes$fbasin <- factor(comm_slopes$Basin.ID)

# Transform EQR changes
pos <- which(comm_slopes$eqr.change > 0)
neg <- which(comm_slopes$eqr.change < 0)
comm_slopes$eqr.sqrt <- sign(comm_slopes$eqr.change) * sqrt(abs(comm_slopes$eqr.change))

# Project WGS84 decimal degrees onto a flat European surface (LCC Europe)
xy_wgs <- comm_slopes[, c("longitude", "latitude")]
xy_lcc <- project(as.matrix(xy_wgs), from = "epsg:4326", to = "epsg:3034")
comm_slopes$lon.lcc <- xy_lcc[, 1]
comm_slopes$lat.lcc <- xy_lcc[, 2]


# Load fitting and predicting functions
source("R Scripts/Trajectories_Analysis functions.R")

### --- Additional metrics --- ###

## Rarefied richness
# Fit
rich.rarefied_gamm <- fit_gamm(
    response = "rich.rarefied",
    data = comm_slopes,
    dist_family = "gaussian"
)
gam.check(rich.rarefied_gamm$gam)
summary(rich.rarefied_gamm$gam)

# Predict
rich.rarefied_rec <- predict_gamm(
    model = rich.rarefied_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0, 0.5, by = 0.1),
    subset_index = pos
)

rich.rarefied_deg <- predict_gamm(
    model = rich.rarefied_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0.5, 1, by = 0.1),
    subset_index = neg
)
rich.rarefied_predicted <- merge(rich.rarefied_rec, rich.rarefied_deg, all = TRUE)


## Coverage-based richness
# Fit
rich.coverage_gamm <- fit_gamm(
    response = "rich.coverage",
    data = comm_slopes,
    dist_family = "gaussian"
)
gam.check(rich.coverage_gamm$gam)
summary(rich.coverage_gamm$gam)

# Predict
rich.coverage_rec <- predict_gamm(
    model = rich.coverage_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0, 0.5, by = 0.1),
    subset_index = pos
)

rich.coverage_deg <- predict_gamm(
    model = rich.coverage_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0.5, 1, by = 0.1),
    subset_index = neg
)
rich.coverage_predicted <- merge(rich.coverage_rec, rich.coverage_deg, all = TRUE)


### --- Export additional metrics predictions --- ###
# Combine predictions
dfs <- list(
    rich.rarefied = rich.rarefied_predicted,
    rich.coverage = rich.coverage_predicted
)
dfs <- purrr::map2(
  dfs,
  names(dfs),
  ~ mutate(.x, response = .y)
)
add_predictions <- bind_rows(dfs)
add_predictions <- add_predictions[, -c(3:8)]
add_predictions <- add_predictions |> relocate(response)

# Export
write.csv(add_predictions, "03_Trajectories_gamm predictions additional.csv", row.names = FALSE)


### --- EQR circularity --- ###

# Define factors
comm_slopes <- comm_slopes |> mutate(
  fric_circularity = as.factor(comm_slopes$richness_circularity),
  fab_circularity = as.factor(comm_slopes$ab_circularity),
  ftaxa_circularity = as.factor(comm_slopes$taxa_circularity),
  ftrait_circularity = as.factor(comm_slopes$trait_circularity)
)

# Number of observations per circularity category
comm_slopes |>
  summarise(
    across(
      c(
        fric_circularity, fab_circularity,
        ftaxa_circularity, ftrait_circularity
      ),
      list(
        "TRUE" = ~ sum(. == TRUE, na.rm = TRUE),
        "FALSE" = ~ sum(. == FALSE, na.rm = TRUE)
      )
    )
  )

# Refit models with factor accounting for circularity

# Richness
# Fit
richness_circular <- fit_gamm_factor(
  response = "richness",
  data = comm_slopes,
  dist_family = "gaussian",
  by_factor = "fric_circularity"
)
gam.check(richness_circular$gam)
summary(richness_circular$gam)

## Re-adjust the range of eqr.sqrt for sites without circularity
richness_false <- comm_slopes[comm_slopes$fric_circularity == FALSE, ]
pos_false <- which(richness_false$eqr.change > 0)
neg_false <- which(richness_false$eqr.change < 0)


# Predict
richness_circular_rec <- predict_gamm_factor(
  model = richness_circular,
  data = richness_false,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos_false,
  by_factor = "fric_circularity"
)

richness_circular_deg <- predict_gamm_factor(
  model = richness_circular,
  data = richness_false,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg_false,
  by_factor = "fric_circularity"
)
richness_circular_predicted <- merge(richness_circular_rec, richness_circular_deg, all = TRUE)


# Abundance
# Fit
ab_circular <- fit_gamm_factor(
  response = "ab",
  data = comm_slopes,
  dist_family = "gaussian",
  by_factor = "fab_circularity"
)
gam.check(ab_circular$gam)
summary(ab_circular$gam)

# Predict

## Re-adjust the range of eqr.sqrt for sites without circularity
ab_false <- comm_slopes[comm_slopes$fab_circularity == FALSE, ]
pos_false <- which(ab_false$eqr.change > 0)
neg_false <- which(ab_false$eqr.change < 0)

ab_circular_rec <- predict_gamm_factor(
  model = ab_circular,
  data = ab_false,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos_false,
  by_factor = "fab_circularity"
)

ab_circular_deg <- predict_gamm_factor(
  model = ab_circular,
  data = ab_false,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg_false,
  by_factor = "fab_circularity"
)
ab_circular_predicted <- merge(ab_circular_rec, ab_circular_deg, all = TRUE)


# Taxon composition
# Fit
taxa_circular <- fit_gamm_factor(
  response = "taxa.comp",
  data = comm_slopes,
  dist_family = "gaussian",
  by_factor = "ftaxa_circularity"
)
gam.check(taxa_circular$gam)
summary(taxa_circular$gam)

# Predict

## Re-adjust the range of eqr.sqrt for sites without circularity
taxa_false <- comm_slopes[comm_slopes$ftaxa_circularity == FALSE, ]
pos_false <- which(taxa_false$eqr.change > 0)
neg_false <- which(taxa_false$eqr.change < 0)

taxa_circular_rec <- predict_gamm_factor(
  model = taxa_circular,
  data = taxa_false,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos_false,
  by_factor = "ftaxa_circularity"
)

taxa_circular_deg <- predict_gamm_factor(
  model = taxa_circular,
  data = taxa_false,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg_false,
  by_factor = "ftaxa_circularity"
)
taxa_circular_predicted <- merge(taxa_circular_rec, taxa_circular_deg, all = TRUE)


# Trait composition
# Fit
trait_circular <- fit_gamm_factor(
  response = "trait.comp",
  data = comm_slopes,
  dist_family = "gaussian",
  by_factor = "ftrait_circularity"
)
gam.check(trait_circular$gam)
summary(trait_circular$gam)

# Predict

## Re-adjust the range of eqr.sqrt for sites without circularity
trait_false <- comm_slopes[comm_slopes$ftrait_circularity == FALSE, ]
pos_false <- which(trait_false$eqr.change > 0)
neg_false <- which(trait_false$eqr.change < 0)

trait_circular_rec <- predict_gamm_factor(
  model = trait_circular,
  data = trait_false,
  eqr_start_vals = seq(0, 0.5, by = 0.1),
  subset_index = pos_false,
  by_factor = "ftrait_circularity"
)

trait_circular_deg <- predict_gamm_factor(
  model = trait_circular,
  data = trait_false,
  eqr_start_vals = seq(0.5, 1, by = 0.1),
  subset_index = neg_false,
  by_factor = "ftrait_circularity"
)
trait_circular_predicted <- merge(trait_circular_rec, trait_circular_deg, all = TRUE)


### --- Export predictions --- ###
# Combine predictions
dfs <- list(
  h0    = h0_circular_predicted,
  ab    = ab_circular_predicted,
  taxa  = taxa_circular_predicted,
  trait = trait_circular_predicted
) |>
  lapply(\(x) { names(x)[3] <- "fcircularity"; x })

dfs <- purrr::map2(
  dfs,
  names(dfs),
  ~ mutate(.x, response = .y)
)
circular_predictions <- bind_rows(dfs)
circular_predictions <- circular_predictions[, -c(4:9)]
circular_predictions <- circular_predictions |> relocate(response)

# Export
write.csv(circular_predictions, "04_Trajectories_gamm predictions circularity.csv", row.names = FALSE)


### --- Taxonomic resolution --- ###

# Define factor
comm_slopes$fid.level <- as.factor(comm_slopes$id.level)

# Repeat main analyses with taxonomic resolution "by" interaction

## Richness
# Fit
richness_id_gamm <- fit_gamm_factor(
    response = "richness",
    data = comm_slopes,
    dist_family = "gaussian"
)
gam.check(richness_id_gamm$gam)
summary(richness_id_gamm$gam)

# Predict
richness_id_rec <- predict_gamm_factor(
    model = richness_id_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0, 0.5, 0.1),
    subset_index = pos,
    by_factor = "fid.level"
)

richness_id_deg <- predict_gamm_factor(
    model = richness_id_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0.5, 1, 0.1),
    subset_index = neg,
    by_factor = "fid.level"
)
richness_id_predicted <- merge(richness_id_rec, richness_id_deg, all = TRUE)


## Taxonomic composition
# Fit
taxa_id_gamm <- fit_gamm_factor(
    response = "taxa.comp",
    data = comm_slopes,
    dist_family = "gaussian"
)
gam.check(taxa_id_gamm$gam)
summary(taxa_id_gamm$gam)

# Predict
taxa_id_rec <- predict_gamm_factor(
    model = taxa_id_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0, 0.5, 0.1),
    subset_index = pos,
    by_factor = "fid.level"
)

taxa_id_deg <- predict_gamm_factor(
    model = taxa_id_gamm,
    data = comm_slopes,
    eqr_start_vals = seq(0.5, 1, 0.1),
    subset_index = neg,
    by_factor = "fid.level"
)
taxa_id_predicted <- merge(taxa_id_rec, taxa_id_deg, all = TRUE)


### --- Export predictions from by factor models --- ###
# Combine predictions
dfs <- list(
  h0 = h0_id_predicted,
  taxa = taxa_id_predicted
)
dfs <- purrr::map2(
  dfs,
  names(dfs),
  ~ mutate(.x, response = .y)
)
factor_predictions <- bind_rows(dfs)
factor_predictions <- factor_predictions[, -c(4:9)]
factor_predictions <- factor_predictions |> relocate(response)

# Export
write.csv(factor_predictions, "05_Trajectories_gamm predictions factor.csv", row.names = FALSE)


### --- Number of sites per country --- ###

# Re sample countries with >= 100 sites to get a more even number of sites
sort(table(comm_slopes$country))
set.seed(123)
n_boot <- 1000

# Identify large countries (>=100 sites)
large_countries <- names(which(table(comm_slopes$country) >= 100))
small_countries <- names(which(table(comm_slopes$country) < 100))

# Create vectors to predict along
eqr_starts_rec <- seq(0, 0.5, by = 0.1)
eqr_starts_deg <- seq(0.5, 1, by = 0.1)

rec_labels <- paste0("pred.r", seq_along(eqr_starts_rec))
deg_labels <- paste0("pred.d", seq_along(eqr_starts_deg))

# Generate eqr.sqrt values for predictions
x_rec <- seq(
    from = min(comm_slopes$eqr.sqrt[pos]),
    to = max(comm_slopes$eqr.sqrt[pos]), length = 1000
)

x_deg <- seq(
    from = min(comm_slopes$eqr.sqrt[neg]),
    to = max(comm_slopes$eqr.sqrt[neg]), length = 1000
)


# - Richness loop
richness_boot_results <- bootstrap_gamm(
    response_var = "richness",
    data = comm_slopes,
    large_countries = large_countries,
    small_countries = small_countries,
    n_boot = n_boot,
    x_rec = x_rec,
    x_deg = x_deg
)

# Store predictions across bootstrap iterations
# Recovery
richness_boot_rec_fit <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)
richness_boot_rec_se <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)

for (i in 1:n_boot) {
  preds <- richness_boot_results$all_preds_rec[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_rec)) {
    richness_boot_rec_fit[, j, i] <- richness_boot_results$all_preds_rec[[i]][[j]]$fit
    richness_boot_rec_se[, j, i] <- richness_boot_results$all_preds_rec[[i]][[j]]$se.fit
  }
}

# Degrading
richness_boot_deg_fit <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)
richness_boot_deg_se <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)

for (i in 1:n_boot) {
  preds <- richness_boot_results$all_preds_deg[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_deg)) {
    richness_boot_deg_fit[, j, i] <- richness_boot_results$all_preds_deg[[i]][[j]]$fit
    richness_boot_deg_se[, j, i] <- richness_boot_results$all_preds_deg[[i]][[j]]$se.fit
  }
}

# Summarize predictions across bootstrap iterations
rec_summary <- data.frame(
    eqr.sqrt = rep(x_rec, times = length(eqr_starts_rec)),
    eqr.start = rep(rec_labels, each = length(x_rec)),
    fit_mean = apply(richness_boot_rec_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(richness_boot_rec_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

deg_summary <- data.frame(
    eqr.sqrt = rep(x_deg, times = length(eqr_starts_deg)),
    eqr.start = rep(deg_labels, each = length(x_deg)),
    fit_mean = apply(richness_boot_deg_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(richness_boot_deg_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

richness_boot <- rbind(rec_summary, deg_summary)
richness_boot <- richness_boot |>
    mutate(
        ci.min = fit_mean - (1.96 * se_mean),
        ci.max = fit_mean + (1.96 * se_mean)
    )

# Average R2adj across bootstrap iterations
richness_r2adj <- richness_boot_results$r2adj
mean(unlist(richness_r2adj, na.rm = TRUE))


# - Abundance loop
ab_boot_results <- bootstrap_gamm(
    response_var = "ab",
    data = comm_slopes,
    large_countries = large_countries,
    small_countries = small_countries,
    n_boot = n_boot,
    x_rec = x_rec,
    x_deg = x_deg
)

# Store predictions across bootstrap iterations
# Recovery
ab_boot_rec_fit <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)
ab_boot_rec_se <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)

for (i in 1:n_boot) {
  preds <- ab_boot_results$all_preds_rec[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_rec)) {
    ab_boot_rec_fit[, j, i] <- ab_boot_results$all_preds_rec[[i]][[j]]$fit
    ab_boot_rec_se[, j, i] <- ab_boot_results$all_preds_rec[[i]][[j]]$se.fit
  }
}

# Degrading
ab_boot_deg_fit <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)
ab_boot_deg_se <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)

for (i in 1:n_boot) {
  preds <- ab_boot_results$all_preds_deg[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_deg)) {
    ab_boot_deg_fit[, j, i] <- ab_boot_results$all_preds_deg[[i]][[j]]$fit
    ab_boot_deg_se[, j, i] <- ab_boot_results$all_preds_deg[[i]][[j]]$se.fit
  }
}

# Summarize predictions across bootstrap iterations to plot
rec_summary <- data.frame(
    eqr.sqrt = rep(x_rec, times = length(eqr_starts_rec)),
    eqr.start = rep(rec_labels, each = length(x_rec)),
    fit_mean = apply(ab_boot_rec_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(ab_boot_rec_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

deg_summary <- data.frame(
    eqr.sqrt = rep(x_deg, times = length(eqr_starts_deg)),
    eqr.start = rep(deg_labels, each = length(x_deg)),
    fit_mean = apply(ab_boot_deg_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(ab_boot_deg_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

ab_boot <- rbind(rec_summary, deg_summary)
ab_boot <- ab_boot |>
    mutate(
        ci.min = fit_mean - (1.96 * se_mean),
        ci.max = fit_mean + (1.96 * se_mean)
    )

# Average R2adj across bootstrap iterations
ab_r2adj <- ab_boot_results$r2adj
mean(unlist(ab_r2adj, na.rm = TRUE))


# - Taxonomic composition loop
taxa_boot_results <- bootstrap_gamm(
    response_var = "taxa.comp",
    data = comm_slopes,
    large_countries = large_countries,
    small_countries = small_countries,
    n_boot = n_boot,
    x_rec = x_rec,
    x_deg = x_deg
)

# Store predictions across bootstrap iterations
# Recovery
taxa_boot_rec_fit <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)
taxa_boot_rec_se <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)

for (i in 1:n_boot) {
  preds <- taxa_boot_results$all_preds_rec[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_rec)) {
    taxa_boot_rec_fit[, j, i] <- taxa_boot_results$all_preds_rec[[i]][[j]]$fit
    taxa_boot_rec_se[, j, i] <- taxa_boot_results$all_preds_rec[[i]][[j]]$se.fit
  }
}

# Degrading
taxa_boot_deg_fit <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)
taxa_boot_deg_se <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)

for (i in 1:n_boot) {
  preds <- taxa_boot_results$all_preds_deg[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_deg)) {
    taxa_boot_deg_fit[, j, i] <- taxa_boot_results$all_preds_deg[[i]][[j]]$fit
    taxa_boot_deg_se[, j, i] <- taxa_boot_results$all_preds_deg[[i]][[j]]$se.fit
  }
}

# Summarize predictions across bootstrap iterations to plot
rec_summary <- data.frame(
    eqr.sqrt = rep(x_rec, times = length(eqr_starts_rec)),
    eqr.start = rep(rec_labels, each = length(x_rec)),
    fit_mean = apply(taxa_boot_rec_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(taxa_boot_rec_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

deg_summary <- data.frame(
    eqr.sqrt = rep(x_deg, times = length(eqr_starts_deg)),
    eqr.start = rep(deg_labels, each = length(x_deg)),
    fit_mean = apply(taxa_boot_deg_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(taxa_boot_deg_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

taxa_boot <- rbind(rec_summary, deg_summary)
taxa_boot <- taxa_boot |>
    mutate(
        ci.min = fit_mean - (1.96 * se_mean),
        ci.max = fit_mean + (1.96 * se_mean)
    )

# Average R2adj across bootstrap iterations
taxa_r2adj <- taxa_boot_results$r2adj
mean(unlist(taxa_r2adj, na.rm = TRUE))


# - Trait composition loop
trait_boot_results <- bootstrap_gamm(
    response_var = "trait.comp",
    data = comm_slopes,
    large_countries = large_countries,
    small_countries = small_countries,
    n_boot = n_boot,
    x_rec = x_rec,
    x_deg = x_deg
)

# Store predictions across bootstrap iterations
# Recovery
trait_boot_rec_fit <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)
trait_boot_rec_se <- array(
    NA,
    dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)

for (i in 1:n_boot) {
  preds <- trait_boot_results$all_preds_rec[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_rec)) {
    trait_boot_rec_fit[, j, i] <- trait_boot_results$all_preds_rec[[i]][[j]]$fit
    trait_boot_rec_se[, j, i] <- trait_boot_results$all_preds_rec[[i]][[j]]$se.fit
  }
}

# Degrading
trait_boot_deg_fit <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)
trait_boot_deg_se <- array(
    NA,
    dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)

for (i in 1:n_boot) {
  preds <- trait_boot_results$all_preds_deg[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_deg)) {
    trait_boot_deg_fit[, j, i] <- trait_boot_results$all_preds_deg[[i]][[j]]$fit
    trait_boot_deg_se[, j, i] <- trait_boot_results$all_preds_deg[[i]][[j]]$se.fit
  }
}

# Summarize predictions across bootstrap iterations to plot
rec_summary <- data.frame(
    eqr.sqrt = rep(x_rec, times = length(eqr_starts_rec)),
    eqr.start = rep(rec_labels, each = length(x_rec)),
    fit_mean = apply(trait_boot_rec_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(trait_boot_rec_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

deg_summary <- data.frame(
    eqr.sqrt = rep(x_deg, times = length(eqr_starts_deg)),
    eqr.start = rep(deg_labels, each = length(x_deg)),
    fit_mean = apply(trait_boot_deg_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
    se_mean = sqrt(apply(trait_boot_deg_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

trait_boot <- rbind(rec_summary, deg_summary)
trait_boot <- trait_boot |>
    mutate(
        ci.min = fit_mean - (1.96 * se_mean),
        ci.max = fit_mean + (1.96 * se_mean)
    )

# Average R2adj across bootstrap iterations
trait_r2adj <- trait_boot_results$r2adj
mean(unlist(trait_r2adj, na.rm = TRUE))


### --- Export bootstrap predictions --- ###
# Combine bootstrap predictions
dfs <- list(
    richness_boot = richness_boot,
    ab_boot = ab_boot,
    taxa_boot = taxa_boot,
    trait_boot = trait_boot
)
dfs <- purrr::map2(
  dfs,
  names(dfs),
  ~ mutate(.x, response = .y)
)
boot_predictions <- bind_rows(dfs)

# Tranform eqr.start labels into values for plotting
names(boot_predictions)[names(boot_predictions) == "eqr.start"] <- "eqr.start.label"
boot_predictions$eqr.start <- c(
  pred.r1 = 0, pred.r2 = 0.1, pred.r3 = 0.2,
  pred.r4 = 0.3, pred.r5 = 0.4, pred.r6 = 0.5,
  pred.d1 = 0.5, pred.d2 = 0.6, pred.d3 = 0.7,
  pred.d4 = 0.8, pred.d5 = 0.9, pred.d6 = 1
)[boot_predictions$eqr.start.label]

boot_predictions <- boot_predictions[, -3]
boot_predictions <- boot_predictions |> relocate(response, eqr.sqrt, eqr.start)

# Export
write.csv(boot_predictions, "06_Trajectories_gamm predictions boot.csv", row.names = FALSE)


### --- Abundance robustness check --- ###

# Re sample recovering sites to get a more even number of sites
range(comm_slopes$eqr.change)
length(comm_slopes$eqr.change[pos]) # 2677
length(comm_slopes$eqr.change[neg]) # 1369
set.seed(123)
n_boot <- 1000

# Split into recovering and degrading sites
pos_sites <- comm_slopes$Unique.ID[pos]
neg_sites <- comm_slopes$Unique.ID[neg]

# Create vectors to predict along
eqr_starts_rec <- seq(0, 0.5, by = 0.1)
eqr_starts_deg <- seq(0.5, 1, by = 0.1)

rec_labels <- paste0("pred.r", seq_along(eqr_starts_rec))
deg_labels <- paste0("pred.d", seq_along(eqr_starts_deg))

# Generate eqr.sqrt values for predictions
x_rec <- seq(
  from = min(comm_slopes$eqr.sqrt[pos]),
  to = max(comm_slopes$eqr.sqrt[pos]), length = 1000
)

x_deg <- seq(
  from = min(comm_slopes$eqr.sqrt[neg]),
  to = max(comm_slopes$eqr.sqrt[neg]), length = 1000
)

# - Abundance loop
ab_resample_results <- resample_gamm(
  response_var = "ab",
  data = comm_slopes,
  pos_sites = pos_sites,
  neg_sites = neg_sites,
  n_boot = n_boot,
  x_rec = x_rec,
  x_deg = x_deg
)

# Store predictions across resampling iterations
# Recovery
ab_resample_rec_fit <- array(
  NA,
  dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)
ab_resample_rec_se <- array(
  NA,
  dim = c(length(x_rec), length(eqr_starts_rec), n_boot)
)

for (i in 1:n_boot) {
  preds <- ab_resample_results$all_preds_rec[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_rec)) {
    ab_resample_rec_fit[, j, i] <- ab_resample_results$all_preds_rec[[i]][[j]]$fit
    ab_resample_rec_se[, j, i] <- ab_resample_results$all_preds_rec[[i]][[j]]$se.fit
  }
}

# Degrading
ab_resample_deg_fit <- array(
  NA,
  dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)
ab_resample_deg_se <- array(
  NA,
  dim = c(length(x_deg), length(eqr_starts_deg), n_boot)
)

for (i in 1:n_boot) {
  preds <- ab_resample_results$all_preds_deg[[i]]
  if (is.null(preds)) {
    next
  }
  for (j in seq_along(eqr_starts_deg)) {
    ab_resample_deg_fit[, j, i] <- ab_resample_results$all_preds_deg[[i]][[j]]$fit
    ab_resample_deg_se[, j, i] <- ab_resample_results$all_preds_deg[[i]][[j]]$se.fit
  }
}

# Summarize predictions across bootstrap iterations to plot
rec_summary <- data.frame(
  eqr.sqrt = rep(x_rec, times = length(eqr_starts_rec)),
  eqr.start = rep(rec_labels, each = length(x_rec)),
  fit_mean = apply(ab_resample_rec_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
  se_mean = sqrt(apply(ab_resample_rec_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

deg_summary <- data.frame(
  eqr.sqrt = rep(x_deg, times = length(eqr_starts_deg)),
  eqr.start = rep(deg_labels, each = length(x_deg)),
  fit_mean = apply(ab_resample_deg_fit, c(1, 2), mean, na.rm = TRUE) |> as.vector(),
  se_mean = sqrt(apply(ab_resample_deg_se^2, c(1, 2), mean, na.rm = TRUE)) |> as.vector()
)

ab_resample <- rbind(rec_summary, deg_summary)
ab_resample <- ab_resample |>
  mutate(
    ci.min = fit_mean - (1.96 * se_mean),
    ci.max = fit_mean + (1.96 * se_mean)
  )

# Average R2adj across bootstrap iterations
ab_resample_r2adj <- ab_resample_results$r2adj
mean(unlist(ab_resample_r2adj, na.rm = TRUE))

# Tranform eqr.start labels into values for plotting
names(ab_resample)[names(ab_resample) == "eqr.start"] <- "eqr.start.label"
ab_resample$eqr.start <- c(
  pred.r1 = 0, pred.r2 = 0.1, pred.r3 = 0.2,
  pred.r4 = 0.3, pred.r5 = 0.4, pred.r6 = 0.5,
  pred.d1 = 0.5, pred.d2 = 0.6, pred.d3 = 0.7,
  pred.d4 = 0.8, pred.d5 = 0.9, pred.d6 = 1
)[ab_resample$eqr.start.label]

ab_resample <- ab_resample[, -2]
ab_resample <- ab_resample |> relocate(eqr.sqrt, eqr.start)

write.csv(ab_resample, "07_Trajectories_gamm predictions boot abundance.csv", row.names = FALSE)
