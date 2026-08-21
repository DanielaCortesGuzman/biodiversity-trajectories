##### ----- Community trajectories analysis functions ----- #####

# This script contains the functions to perform the analysis
# of community trajectories in Trajectories_Main analysis.R
# and Trajectories_Supplementary analysis.R


# Define function for GAMM fitting
fit_gamm <- function(response,
                     data,
                     dist_family,
                     add_fixed = NULL) {
    # Based terms
    base_terms <- c(
        "te(eqr.sqrt, eqr.start)",
        "s(fcountry, bs = 're')",
        "s(fbasin, bs = 're')",
        "lon.lcc", "lat.lcc", "dur.yr", "st.yr"
    )

    # Additional fixed terms if provided
    if (!is.null(add_fixed)) {
        all_terms <- c(base_terms, add_fixed)
    } else {
        all_terms <- base_terms
    }

    # gamm formula
    form <- as.formula(
        paste(
            response, "~", paste(all_terms, collapse = " + ")
        )
    )

    # model
    gamm(
        formula = form,
        data = data,
        family = dist_family,
        method = "REML"
    )
}


# Define prediction function for GAMM models
predict_gamm <- function(model,
                         data,
                         eqr_start_vals,
                         subset_index,
                         n_points = 1000,
                         exclude_term = c("s(fcountry)", "s(fbasin)"),
                         add_fixed = NULL) {
    # Eqr.sqrt sequence
    x_seq <- seq(
        from = min(data$eqr.sqrt[subset_index]),
        to = max(data$eqr.sqrt[subset_index]),
        length.out = n_points
    )

    # Prediction dataframe
    base_df <- expand.grid(
        eqr.sqrt = x_seq,
        eqr.start = eqr_start_vals
    )

    # Add fixed covariates
    newdata <- base_df |>
        mutate(
            fcountry = data$fcountry[1],
            fbasin = data$fbasin[1],
            lat.lcc = mean(data$lat.lcc),
            lon.lcc = mean(data$lon.lcc),
            dur.yr = mean(data$dur.yr),
            st.yr = mean(data$st.yr)
        )

    # Add additional fixed predictors, when needed
    if (!is.null(add_fixed)) {
        for (var in add_fixed) {
            if (is.numeric(data[[var]])) {
                newdata[[var]] <- mean(data[[var]], na.rm = TRUE)
            } else if (is.factor(data[[var]]) || is.character(data[[var]])) {
                newdata[[var]] <- data[[var]][1]
            } else {
                stop(paste("Unsopported variable type for: ", var))
            }
        }
    }

    # Predict
    preds <- predict(
        model$gam,
        newdata = newdata,
        type = "response",
        exclude = exclude_term,
        se.fit = TRUE
    )

    # Bind predictions
    result <- newdata |>
        mutate(
            fit    = preds$fit,
            se.fit = preds$se.fit,
            ci.min = fit - 1.96 * se.fit,
            ci.max = fit + 1.96 * se.fit
        )

    return(result)
}


# Define GAMM function without fixed terms (for checking purposes)
fit_gamm_fixed <- function(response,
                           data,
                           dist_family) {
    # gamm form
    form <- as.formula(
        paste(
            response, "~ te(eqr.sqrt, eqr.start) +",
            "s(fcountry, bs = 're') +",
            "s(fbasin, bs = 're')"
        )
    )

    # model
    gamm(
        formula = form,
        data = data,
        family = dist_family,
        method = "REML"
    )
}


# Define function for GAMM with a "by" factor
fit_gamm_factor <- function(response,
                            data,
                            dist_family,
                            by_factor) {
    # gamm form
    form <- as.formula(
        paste(
            response, "~ te(eqr.sqrt, eqr.start, by =", by_factor, ") +",
            "s(fcountry, bs = 're') +",
            "s(fbasin, bs = 're') +",
            "lon.lcc + lat.lcc + dur.yr + st.yr +", by_factor
        )
    )

    # model
    gamm(
        formula = form,
        data = data,
        family = dist_family,
        method = "REML"
    )
}


# Define prediction function for GAMM models with a "by" factor
predict_gamm_factor <- function(model,
                                data,
                                eqr_start_vals,
                                subset_index,
                                by_factor,
                                n_points = 1000,
                                exclude_term = c("s(fcountry)", "s(fbasin)")) {
    # Eqr.sqrt sequence
    x_seq <- seq(
        from = min(data$eqr.sqrt[subset_index]),
        to = max(data$eqr.sqrt[subset_index]),
        length.out = n_points
    )

    # Taxa resolution levels
    factor_levels <- levels(data[[by_factor]])

    # Prediction dataframe
    newdata <- expand.grid(
        eqr.sqrt = x_seq,
        eqr.start = eqr_start_vals,
        level = factor_levels
    )
    names(newdata)[names(newdata) == "level"] <- by_factor
    newdata[[by_factor]] <- factor(
        newdata[[by_factor]],
        levels = factor_levels
    )

    # Add fixed covariates
    newdata <- newdata |>
        mutate(
            fcountry = data$fcountry[1],
            fbasin = data$fbasin[1],
            lat.lcc = mean(data$lat.lcc),
            lon.lcc = mean(data$lon.lcc),
            dur.yr = mean(data$dur.yr),
            st.yr = mean(data$st.yr)
        )

    # Predict
    preds <- predict(
        model$gam,
        newdata = newdata,
        type = "response",
        exclude = exclude_term,
        se.fit = TRUE
    )

    # Bind predictions
    result <- newdata |>
        mutate(
            fit    = preds$fit,
            se.fit = preds$se.fit,
            ci.min = fit - 1.96 * se.fit,
            ci.max = fit + 1.96 * se.fit
        )

    return(result)
}


# Define function to bootstrap GAMM models with the number of sites
bootstrap_gamm <- function(response_var,
                           data,
                           large_countries,
                           small_countries,
                           n_boot,
                           x_rec,
                           x_deg) {
    # lists to store results
    r2adj <- vector("list", n_boot)
    all_preds_rec <- vector("list", n_boot)
    all_preds_deg <- vector("list", n_boot)

    # build formula
    form <- as.formula(
        paste0(
            response_var,
            " ~ te(eqr.sqrt, eqr.start) + ",
            "s(fcountry, bs = 're') +
            s(fbasin, bs = 're') +
            lon.lcc + lat.lcc + dur.yr + st.yr"
        )
    )

    # - loop
    for (i in 1:n_boot) {
        # Resample large countries
        resampled_indices <- c()
        for (ctry in large_countries) {
            id <- which(data$country == ctry)
            resampled_indices <- c(
                resampled_indices,
                sample(id, size = 100, replace = FALSE)
            )
        }

        # Include all small countries
        small_indices <- unlist(
            lapply(
                small_countries,
                function(ctry) which(data$country == ctry)
            )
        )

        # Combine large (resampled) and small
        boot_data <- data[c(resampled_indices, small_indices), ]

        # Fit GAMM
        fit <- tryCatch(
            gamm(
                form,
                data = boot_data,
                family = "gaussian",
                method = "REML"
            ),
            error = function(e) NULL
        )

        if (is.null(fit)) {
            cat("Bootstrap iteration:", i, "failed to fit the model.\n")
            all_preds_rec[[i]] <- NULL
            all_preds_deg[[i]] <- NULL
            r2adj[[i]] <- NA_real_

            next
        }

        # Store adjusted R-squared
        r2adj[[i]] <- summary(fit$gam)$r.sq

        # Recovery predictions
        eqr_starts_rec <- seq(0, 0.5, by = 0.1)
        preds_rec <- list()
        for (j in seq_along(eqr_starts_rec)) {
            preds_rec[[j]] <- predict(fit$gam,
                newdata = data.frame(
                    eqr.sqrt = x_rec,
                    eqr.start = eqr_starts_rec[j],
                    fcountry = boot_data$fcountry[1],
                    fbasin = boot_data$fbasin[1],
                    lat.lcc = mean(boot_data$lat.lcc),
                    lon.lcc = mean(boot_data$lon.lcc),
                    dur.yr = mean(boot_data$dur.yr),
                    st.yr = mean(boot_data$st.yr)
                ),
                type = "response",
                exclude = c("s(fcountry)", "s(fbasin)"),
                se.fit = TRUE
            )
        }
        all_preds_rec[[i]] <- preds_rec

        # Degrading predictions
        eqr_starts_deg <- seq(0.5, 1, by = 0.1)
        preds_deg <- list()
        for (j in seq_along(eqr_starts_deg)) {
            preds_deg[[j]] <- predict(fit$gam,
                newdata = data.frame(
                    eqr.sqrt = x_deg,
                    eqr.start = eqr_starts_deg[j],
                    fcountry = boot_data$fcountry[1],
                    fbasin = boot_data$fbasin[1],
                    lat.lcc = mean(boot_data$lat.lcc),
                    lon.lcc = mean(boot_data$lon.lcc),
                    dur.yr = mean(boot_data$dur.yr),
                    st.yr = mean(boot_data$st.yr)
                ),
                type = "response",
                exclude = c("s(fcountry)", "s(fbasin)"),
                se.fit = TRUE
            )
        }
        all_preds_deg[[i]] <- preds_deg

        cat("Bootstrap iteration:", i, "\n")
    }

    return(list(
        r2adj = r2adj,
        all_preds_rec = all_preds_rec,
        all_preds_deg = all_preds_deg
    ))
}


# Define function to bootstrap GAMM models with even number of degrading and recovering sites
resample_gamm <- function(response_var,
                          data,
                          neg_sites,
                          pos_sites,
                          n_boot,
                          x_rec,
                          x_deg) {
    # lists to store results
    r2adj <- vector("list", n_boot)
    all_preds_rec <- vector("list", n_boot)
    all_preds_deg <- vector("list", n_boot)

    # build formula
    form <- as.formula(
        paste0(
            response_var,
            " ~ te(eqr.sqrt, eqr.start) + ",
            "s(fcountry, bs = 're') +
            s(fbasin, bs = 're') +
            lon.lcc + lat.lcc + dur.yr + st.yr"
        )
    )

    # - loop
    for (i in 1:n_boot) {
        # Resample recovering sites
        id <- which(data$Unique.ID %in% pos_sites)
        resampled_indices <- sample(id, size = 1369, replace = FALSE)

        # Include all degrading sites
        neg_indices <- which(data$Unique.ID %in% neg_sites)

        # Combine recovergin (resampled) and degrading
        boot_data <- data[c(resampled_indices, neg_indices), ]

        # Fit GAMM
        fit <- tryCatch(
            gamm(
                form,
                data = boot_data,
                family = "gaussian",
                method = "REML"
            ),
            error = function(e) NULL
        )

        if (is.null(fit)) {
            cat("Bootstrap iteration:", i, "failed to fit the model.\n")
            all_preds_rec[[i]] <- NULL
            all_preds_deg[[i]] <- NULL
            r2adj[[i]] <- NA_real_
            next
        }

        # Store adjusted R-squared
        r2adj[[i]] <- summary(fit$gam)$r.sq

        # Recovery predictions
        eqr_starts_rec <- seq(0, 0.5, by = 0.1)
        preds_rec <- list()
        for (j in seq_along(eqr_starts_rec)) {
            preds_rec[[j]] <- predict(fit$gam,
                newdata = data.frame(
                    eqr.sqrt = x_rec,
                    eqr.start = eqr_starts_rec[j],
                    fcountry = boot_data$fcountry[1],
                    fbasin = boot_data$fbasin[1],
                    lat.lcc = mean(boot_data$lat.lcc),
                    lon.lcc = mean(boot_data$lon.lcc),
                    dur.yr = mean(boot_data$dur.yr),
                    st.yr = mean(boot_data$st.yr)
                ),
                type = "response",
                exclude = c("s(fcountry)", "s(fbasin)"),
                se.fit = TRUE
            )
        }
        all_preds_rec[[i]] <- preds_rec

        # Degrading predictions
        eqr_starts_deg <- seq(0.5, 1, by = 0.1)
        preds_deg <- list()
        for (j in seq_along(eqr_starts_deg)) {
            preds_deg[[j]] <- predict(fit$gam,
                newdata = data.frame(
                    eqr.sqrt = x_deg,
                    eqr.start = eqr_starts_deg[j],
                    fcountry = boot_data$fcountry[1],
                    fbasin = boot_data$fbasin[1],
                    lat.lcc = mean(boot_data$lat.lcc),
                    lon.lcc = mean(boot_data$lon.lcc),
                    dur.yr = mean(boot_data$dur.yr),
                    st.yr = mean(boot_data$st.yr)
                ),
                type = "response",
                exclude = c("s(fcountry)", "s(fbasin)"),
                se.fit = TRUE
            )
        }
        all_preds_deg[[i]] <- preds_deg

        cat("Bootstrap iteration:", i, "\n")
    }

    return(list(
        r2adj = r2adj,
        all_preds_rec = all_preds_rec,
        all_preds_deg = all_preds_deg
    ))
}
