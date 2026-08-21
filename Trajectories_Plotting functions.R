##### ----- Community trajectories plotting functions ----- #####

# This script contains the functions to produce
# visualizations of community trajectories in Trajectories_Plots.R

# Define plotting function with default values
plot_trajectory_trend <- function(df,
                                  fit_var = "fit",
                                  ci_min_var = "ci.min",
                                  ci_max_var = "ci.max",
                                  scale_y = percent_format(scale = 1, accuracy = 1),
                                  y_breaks = c(-8, -4, 0, 4, 8),
                                  y_limits = c(-8, 8),
                                  y_label = expression("Trend (" * year^{
                                      -1
                                  } * ")"),
                                  x_label = "Total change in ecological quality",
                                  plot_title = NULL) {
    ggplot(df) +
        geom_ribbon(
            aes(
                x = eqr.sqrt,
                ymin = .data[[ci_min_var]],
                ymax = .data[[ci_max_var]],
                fill = feqr.start2
            ),
            alpha = 0.1, linetype = 0
        ) +
        geom_line(
            aes(
                x = eqr.sqrt,
                y = .data[[fit_var]],
                color = feqr.start2
            ),
            linewidth = 1
        ) +
        scale_colour_manual(values = eqr_colors) +
        scale_fill_manual(values = eqr_colors) +
        geom_hline(
            yintercept = 0, color = "black",
            linetype = "dashed", linewidth = 1
        ) +
        scale_y_continuous(
            labels = scale_y,
            breaks = y_breaks,
            limits = y_limits
        ) +
        scale_x_continuous(
            breaks = x_breaks, labels = x_labels,
            minor_breaks = minor_breaks,
            limits = x_limits,
            guide = guide_axis(minor.ticks = TRUE)
        ) +
        labs(y = y_label, x = x_label, title = plot_title) +
        theme_bw(base_size = 20) +
        theme(
            text = element_text(size = 20),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks = element_line(linewidth = 0.8, color = "black"),
            axis.ticks.length = unit(0.3, "cm"),
            legend.position = "none",
            axis.minor.ticks.length = rel(1),
            axis.minor.ticks.x.bottom = element_line(colour = "gray")
        )
}


# Define plotting with data points
plot_trajectory_points <- function(df,
                                   var = "var",
                                   fit_var = "fit",
                                   ci_min_var = "ci.min",
                                   ci_max_var = "ci.max",
                                   scale_y = percent_format(scale = 1, accuracy = 1),
                                   y_breaks = c(-8, -4, 0, 4, 8),
                                   y_limits = c(-8, 8),
                                   y_label = expression("Trend (" * year^{
                                       -1
                                   } * ")"),
                                   x_label = "Total change in ecological quality",
                                   plot_title = NULL) {
    ggplot(df) +
        geom_point(
            aes(x = eqr.sqrt, y = .data[[var]]),
            data = comm_slopes,
            size = 1, color = "lightgrey", shape = 20
        ) +
        geom_ribbon(
            aes(
                x = eqr.sqrt,
                ymin = .data[[ci_min_var]],
                ymax = .data[[ci_max_var]],
                fill = feqr.start2
            ),
            alpha = 0.1, linetype = 0
        ) +
        geom_line(
            aes(
                x = eqr.sqrt,
                y = .data[[fit_var]],
                color = feqr.start2
            ),
            linewidth = 1
        ) +
        scale_colour_manual(values = eqr_colors) +
        scale_fill_manual(values = eqr_colors) +
        geom_hline(
            yintercept = 0, color = "black",
            linetype = "dashed", linewidth = 1
        ) +
        scale_y_continuous(
            labels = scale_y,
            breaks = y_breaks,
            limits = y_limits
        ) +
        scale_x_continuous(
            breaks = x_breaks_ext, labels = x_labels_ext,
            minor_breaks = minor_breaks_ext,
            limits = x_limits_ext,
            guide = guide_axis(minor.ticks = TRUE)
        ) +
        labs(y = y_label, x = x_label, title = plot_title) +
        theme_bw(base_size = 20) +
        theme(
            text = element_text(size = 20),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks = element_line(linewidth = 0.8, color = "black"),
            axis.ticks.length = unit(0.3, "cm"),
            legend.position = "none",
            axis.minor.ticks.length = rel(1),
            axis.minor.ticks.x.bottom = element_line(colour = "gray")
        )
}
