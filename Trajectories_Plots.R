###### ----- Community trajectories plots ----- #####

# This script produces visualizations of community trajectories predicted values
# in relationship to EQR changes,
# based on the GAMM analyses performed in Trajectories_Main analysis.R
# and supplementary plots from Trajectories_Supplementary analysis.R

library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(purrr)
library(scico)


# Load plotting functions
source("R Scripts/Trajectories_Plotting functions.R")


# Load prediction files
# List of files
names <- c("main", "add", "circular", "factor", "boot_country", "boot_ab")
files <- lapply(1:length(names), function(x) file.choose())

# Load all files and define factors for plotting
datasets <- lapply(files, function(f) {
    df <- read.csv(f, header = TRUE)

    if ("eqr.start" %in% names(df)) {
        df$feqr.start <- factor(df$eqr.start, levels = sort(unique(df$eqr.start)))
    }

    if ("fid.level" %in% names(df)) {
        df$fid.level <- as.factor(df$fid.level)
    }

    df
})

# Create factor for EQR start values (to distinguish 0.5 between rec and deg)
datasets <- lapply(datasets, function(df) {
    if ("eqr.start" %in% names(df)) {
        df$feqr.start2 <- ifelse(df$eqr.start == 0.5 & df$eqr.sqrt > 0, 0.45, df$eqr.start)
        df$feqr.start2 <- factor(df$feqr.start2, levels = sort(unique(df$feqr.start2)))
    }
    df
})

# Separate individual dataframes
names(datasets) <- names
main <- datasets$main
add <- datasets$add
eqr <- datasets$eqr
factor <- datasets$factor
boot_country <- datasets$boot_country
boot_ab <- datasets$boot_ab

# Define plot features
x_breaks <- c(-0.71, -0.55, -0.39, -0.22, 0, 0.22, 0.39, 0.55, 0.71)
x_labels <- c(
    "-0.5", "-0.3", "-0.15", "-0.05",
    "0", "0.05", "0.15", "0.3", "0.5"
)
x_limits <- c(-0.71, 0.71)
eqr_colors <- c(
    "0" = "#440154",
    "0.1" = "#482878",
    "0.2" = "#31688E",
    "0.3" = "#1F9E89",
    "0.4" = "#6DCD59",
    "0.45" = "#FDE725",
    "0.5" = "#440154",
    "0.6" = "#482878",
    "0.7" = "#31688E",
    "0.8" = "#1F9E89",
    "0.9" = "#6DCD59",
    "1" = "#FDE725"
)

minor_breaks <- c(
    -0.67, -0.63, -0.59, -0.5, -0.45, -0.32,
    0.32, 0.45, 0.5, 0.59, 0.63, 0.67
)


## -- Richness -- ##
# Plot
plot_trajectory_trend(
    df = subset(main, response == "richness"),
    fit_var = "fit",
    ci_min_var = "ci.min",
    ci_max_var = "ci.max",
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-6, -4, -2, 0, 2, 4, 6, 8),
    y_limits = c(-6, 8),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("richness.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Abundance -- ##
plot_trajectory_trend(
    df = subset(main, response == "ab"),
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-4, -2, 0, 2, 4, 6),
    y_limits = c(-5, 6),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("ab.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Functional diversity -- ##
plot_trajectory_trend(
    df = subset(main, response == "func.div"),
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-8, -4, 0, 4, 8, 12),
    y_limits = c(-10, 13),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("func.div.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Trait SES -- ##
plot_trajectory_trend(
    df = subset(main, response == "ses"),
    scale_y = label_number(accuracy = 0.01),
    y_breaks = c(-0.08, -0.04, 0, 0.04, 0.08, 0.12),
    y_limits = c(-0.08, 0.12),
    y_label = expression("Trend (SD " * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("ses.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Taxonomic composition -- ##
range(main$taxa_fit)
plot_trajectory_trend(
    df = subset(main, response == "taxa"),
    scale_y = label_number(accuracy = 0.01),
    y_breaks = c(0.04, 0.06, 0.08, 0.10, 0.12),
    y_limits = c(0.03, 0.13),
    y_label = expression(
        "Compositional change (SD " ~ year^
            {
                -1
            } * ")"
    ),
    x_label = "Total change in ecological quality"
)
ggsave("taxa_comp.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Trait composition -- ##
range(main$trait_fit)
plot_trajectory_trend(
    df = subset(main, response == "trait"),
    scale_y = label_number(accuracy = 0.01),
    y_breaks = c(0.01, 0.02, 0.03, 0.04),
    y_limits = c(0.004, 0.04),
    y_label = expression(
        "Compositional change (SD " ~ year^
            {
                -1
            } * ")"
    ),
    x_label = "Total change in ecological quality"
)
ggsave("trait_comp.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)



### --- Supplementary plots -- ###

## -- Main plots with data points

# Load data
comm_slopes <- read.csv(file.choose(), header = TRUE) # 01_Trajectories_biodiversity slopes.csv

# Extend x-axis limits and breaks for full visualization of data points
comm_slopes$eqr.sqrt <- sign(comm_slopes$eqr.change) * sqrt(abs(comm_slopes$eqr.change))
pos <- which(comm_slopes$eqr.change > 0)
neg <- which(comm_slopes$eqr.change < 0)
range(comm_slopes$eqr.sqrt[pos])
range(comm_slopes$eqr.sqrt[neg])
range(comm_slopes$eqr.change)

x_breaks_ext <- c(-0.71, -0.55, -0.39, -0.22, 0, 0.22, 0.39, 0.55, 0.71, 0.89)
x_labels_ext <- c(
  "-0.5", "-0.3", "-0.15", "-0.05",
  "0", "0.05", "0.15", "0.3", "0.5", "0.8"
)
x_limits_ext <- c(-0.81, 0.91)
minor_breaks_ext <- c(
  -0.77, -0.67, -0.63, -0.59, -0.5, -0.45, -0.32,
  0.32, 0.45, 0.5, 0.59, 0.63, 0.67, 0.77, 0.84
)


## -- Richness -- ##
range(comm_slopes$richness)
plot_trajectory_points(
  df = subset(main, response == "richness"),
  var = "richness",
  scale_y = percent_format(scale = 1, accuracy = 1),
  y_breaks = c(-12, -6, 0, 6, 12, 18),
  y_limits = c(-12, 20),
  y_label = expression("Trend (" * year^{
    -1
  } * ")"),
  x_label = "Total change in ecological quality"
)
ggsave("richness_points.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Abundance -- ##
range(comm_slopes$ab)
plot_trajectory_points(
  df = subset(main, response == "ab"),
  var = "ab",
  scale_y = percent_format(scale = 1, accuracy = 1),
  y_breaks = c(-40, -20, 0, 20, 40, 60, 80),
  y_limits = c(-40, 95),
  y_label = expression("Trend (" * year^{
    -1
  } * ")"),
  x_label = "Total change in ecological quality"
)
ggsave("ab_points.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Functional diversity -- ##
range(comm_slopes$func.div, na.rm = TRUE)
plot_trajectory_points(
  df = subset(main, response == "fd"),
  var = "fd",
  scale_y = percent_format(scale = 1, accuracy = 1),
  y_breaks = c(-20, -10, 0, 10, 20, 30, 40),
  y_limits = c(-23, 40),
  y_label = expression("Trend (" * year^{
    -1
  } * ")"),
  x_label = "Total change in ecological quality"
)
ggsave("fd_points.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Trait SES -- ##
range(comm_slopes$trait.ses, na.rm = TRUE)
plot_trajectory_points(
  df = subset(main, response == "ses"),
  var = "trait.ses",
  scale_y = label_number(accuracy = 0.01),
  y_breaks = c(-0.20, -0.10, 0, 0.10, 0.20, 0.30),
  y_limits = c(-0.26, 0.36),
  y_label = expression("Trend (SD " * year^{
    -1
  } * ")"),
  x_label = "Total change in ecological quality"
)
ggsave("ses_points.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Taxonomic composition -- ##
range(comm_slopes$taxa.comp, na.rm = TRUE)
plot_trajectory_points(
  df = subset(main, response == "taxa"),
  var = "taxa.comp",
  scale_y = label_number(accuracy = 0.01),
  y_breaks = c(0, 0.1, 0.2, 0.3, 0.4),
  y_limits = c(0.00001, 0.4),
  y_label = expression(
    "Compositional change (SD " ~ year^
      {
        -1
      } * ")"
  ),
  x_label = "Total change in ecological quality"
)
ggsave("taxa_comp_points.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Trait composition -- ##
range(comm_slopes$trait.comp, na.rm = TRUE)
plot_trajectory_points(
  df = subset(main, response == "trait"),
  var = "trait.comp",
  scale_y = label_number(accuracy = 0.01),
  y_breaks = c(0, 0.03, 0.06, 0.09, 0.12),
  y_limits = c(0.00001, 0.135),
  y_label = expression(
    "Compositional change (SD " ~ year^
      {
        -1
      } * ")"
  ),
  x_label = "Total change in ecological quality"
)
ggsave("trait_comp_points.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Additional metrics

# - Rarefied richness
plot_trajectory_trend(
    df = subset(add, response == "rich.rarefied"),
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-4, -2, 0, 2, 4, 6, 8),
    y_limits = c(-7, 8),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("rich.rarefied.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


# - Coverage-based richness
plot_trajectory_trend(
    df = subset(add, response == "rich.coverage"),
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8, 10),
    y_limits = c(-10, 10),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("rich.coverage.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Circularity check

# Split into separate dataframes for each level
circular_false <- circular |> filter(fcircularity == FALSE)
circular_true <- circular |> filter(fcircularity == TRUE)

comm_slopes <- comm_slopes |> mutate(
  richness_removed = ifelse(richness_circularity == TRUE, TRUE, FALSE),
  ab_removed = ifelse(ab_circularity == TRUE, TRUE, FALSE),
  taxa_removed = ifelse(taxa_circularity == TRUE, TRUE, FALSE),
  trait_removed = ifelse(trait_circularity == TRUE, TRUE, FALSE)
)

## -- Richness -- ##
# Plot
plot_trajectory_trend(
  df = subset(circular_false, response == "richness"),
  fit_var = "fit",
  ci_min_var = "ci.min",
  ci_max_var = "ci.max",
  scale_y = percent_format(scale = 1, accuracy = 1),
  y_breaks = c(-6, -4, -2, 0, 2, 4, 6, 8),
  y_limits = c(-6, 8),
  y_label = expression("Trend (" * year^{
    -1
  } * ")"),
  x_label = "Total change in ecological quality"
)
ggsave("richness_circularity.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Abundance -- ##
plot_trajectory_trend(
  df = subset(circular_false, response == "ab"),
  fit_var = "fit",
  ci_min_var = "ci.min",
  ci_max_var = "ci.max",
  scale_y = percent_format(scale = 1, accuracy = 1),
  y_breaks = c(-4, -2, 0, 2, 4, 6),
  y_limits = c(-5, 6),
  y_label = expression("Trend (" * year^{
    -1
  } * ")"),
  x_label = "Total change in ecological quality"
)
ggsave("ab_circularity.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Taxonomic composition -- ##
plot_trajectory_trend(
  df = subset(circular_false, response == "taxa"),
  fit_var = "fit",
  ci_min_var = "ci.min",
  ci_max_var = "ci.max",
  scale_y = label_number(accuracy = 0.01),
  y_breaks = c(0.04, 0.06, 0.08, 0.10, 0.12),
  y_limits = c(0.025, 0.13),
  y_label = expression(
    "Compositional change (SD " ~ year^
      {
        -1
      } * ")"
  ),
  x_label = "Total change in ecological quality"
)
ggsave("taxa_circularity.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Trait composition -- ##
plot_trajectory_trend(
  df = subset(circular_false, response == "trait"),
  fit_var = "fit",
  ci_min_var = "ci.min",
  ci_max_var = "ci.max",
  scale_y = label_number(accuracy = 0.01),
  y_breaks = c(0.01, 0.02, 0.03, 0.04),
  y_limits = c(0.003, 0.04),
  y_label = expression(
    "Compositional change (SD " ~ year^
      {
        -1
      } * ")"
  ),
  x_label = "Total change in ecological quality"
)
ggsave("trait_circularity.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)



## -- Taxonomic resolution

# Split into separate dataframes for each taxonomic resolution
factor_mixed <- factor |> filter(fid.level == "Mixed")
factor_family <- factor |> filter(fid.level == "Family")

# - Richness
plot_trajectory_trend(
    df = subset(factor_mixed, response == "richness"),
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-8, -4, 0, 4, 8),
    y_limits = c(-8, 11),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("richness_mixed.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)

plot_trajectory_trend(
    df = subset(factor_family, response == "richness"),
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-8, -4, 0, 4, 8),
    y_limits = c(-8, 11),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("richness_family.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


# - Taxonomic composition
plot_trajectory_trend(
    df = subset(factor_mixed, response == "taxa"),
    scale_y = label_number(accuracy = 0.01),
    y_breaks = c(0.04, 0.06, 0.08, 0.10, 0.12),
    y_limits = c(0.02, 0.13),
    y_label = expression(
        "Compositional change (SD " ~ year^
            {
                -1
            } * ")"
    ),
    x_label = "Total change in ecological quality"
)
ggsave("taxa_mixed.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)

plot_trajectory_trend(
    df = subset(factor_family, response == "taxa"),
    scale_y = label_number(accuracy = 0.01),
    y_breaks = c(0.04, 0.06, 0.08, 0.10, 0.12),
    y_limits = c(0.02, 0.13),
    y_label = expression(
        "Compositional change (SD " ~ year^
            {
                -1
            } * ")"
    ),
    x_label = "Total change in ecological quality"
)
ggsave("taxa_family.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Number of sites per country

# - Richness
plot_trajectory_trend(
    df = subset(boot_country, response == "richness"),
    fit_var = "fit_mean",
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-6, -4, -2, 0, 2, 4, 6, 8),
    y_limits = c(-6, 8),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("richness_boot.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


# - Abundance
plot_trajectory_trend(
    df = subset(boot_country, response == "ab"),
    fit_var = "fit_mean",
    scale_y = percent_format(scale = 1, accuracy = 1),
    y_breaks = c(-4, -2, 0, 2, 4, 6),
    y_limits = c(-5, 6),
    y_label = expression("Trend (" * year^{
        -1
    } * ")"),
    x_label = "Total change in ecological quality"
)
ggsave("ab_boot.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


# - Taxonomic composition
plot_trajectory_trend(
    df = subset(boot_country, response == "taxa"),
    fit_var = "fit_mean",scale_y = label_number(accuracy = 0.01),
    y_breaks = c(0.04, 0.06, 0.08, 0.10, 0.12),
    y_limits = c(0.03, 0.13),
    y_label = expression(
        "Compositional change (SD " ~ year^
            {
                -1
            } * ")"
    ),
    x_label = "Total change in ecological quality"
)
ggsave("taxa_boot.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


# - Trait composition
plot_trajectory_trend(
    df = subset(boot_country, response == "trait"),
    fit_var = "fit_mean",
    scale_y = label_number(accuracy = 0.01),
    y_breaks = c(0.01, 0.02, 0.03, 0.04),
    y_limits = c(0.004, 0.04),
    y_label = expression(
        "Compositional change (SD " ~ year^
            {
                -1
            } * ")"
    ),
    x_label = "Total change in ecological quality"
)
ggsave("trait_boot.svg",
    device = svglite::svglite, width = 170, height = 120, units = "mm"
)


## -- Resampling check: degrading and recovering sites

# - Abundance
plot_trajectory_trend(
  df = boot_ab,
  fit_var = "fit_mean",
  ci_min_var = "ci.min",
  ci_max_var = "ci.max",
  scale_y = percent_format(scale = 1, accuracy = 1),
  y_breaks = c(-4, -2, 0, 2, 4, 6),
  y_limits = c(-5, 6),
  y_label = expression("Trend (" * year^{
    -1
  } * ")"),
  x_label = "Total change in ecological quality"
)
ggsave("ab_resampled.svg",
       device = svglite::svglite, width = 170, height = 120, units = "mm"
)