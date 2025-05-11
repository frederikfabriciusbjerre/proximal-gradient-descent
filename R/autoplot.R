library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

autoplot.pcd_out <- function(x, tol = NULL) {
  if (is.null(tol)) {
    warning("`tol` is NULL, using default value of 1e-15")
    tol <- 1e-15
  }

  # pivot to long
  df <- as_tibble(x$path) |>
    mutate(iter = seq_len(n())) |>
    pivot_longer(
      cols = -iter,
      names_to = "term",
      values_to = "estimate"
    )

  # get each term's final value
  final_vals <- df |>
    group_by(term) |>
    summarize(final = last(estimate), .groups = "drop")

  # find first iteration where estimate is within tol of final
  conv_pts <- df |>
    inner_join(final_vals, by = "term") |>
    group_by(term) |>
    filter(abs(estimate - final) < tol) |>
    slice_head(n = 1) |>
    ungroup()

  # plot lines + circles at convergence
  ggplot(df, aes(x = iter, y = estimate, colour = term)) +
    geom_line() +
    geom_point(
      data = conv_pts,
      aes(x = iter, y = estimate, colour = term),
      shape = 21, # hollow circle
      size = 3,
      stroke = 1.2
    ) +
    labs(
      x = "Iteration",
      y = expression(hat(beta)),
      colour = "Coefficient"
    ) +
    theme_minimal()
}
