#' Plot convergence plots for pcd objects
#'
#' @param x A `pcd` object.
#' @param tol Tolerance for convergence. Default is `NULL`,
#' which will use a default value of `1e-12`.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom dplyr mutate filter group_by summarize inner_join slice_head ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @export
autoplot.pcd_out <- function(x, tol = NULL, legend = TRUE) {
  if (is.null(tol)) {
    warning("`tol` is NULL, using default value of 1e-12")
    tol <- 1e-12
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
  # now build the plot in four layers:
  p <- ggplot() +
    # all non-intercept lines, coloured by term
    geom_line(
      data = df |> filter(term != "(Intercept)"),
      aes(x = iter, y = estimate, colour = term)
    ) +
    # intercept line, black & dashed
    geom_line(
      data = df |> filter(term == "(Intercept)"),
      aes(x = iter, y = estimate),
      colour = "black",
      linetype = "dashed",
      size = 1
    ) +
    # non-intercept convergence points
    geom_point(
      data = conv_pts |> filter(term != "(Intercept)"),
      aes(x = iter, y = estimate, colour = term),
      shape = 21, size = 3, stroke = 1.2
    ) +
    # intercept convergence point, black outline
    geom_point(
      data = conv_pts |> filter(term == "(Intercept)"),
      aes(x = iter, y = estimate),
      shape = 21, size = 3, stroke = 1.2,
      colour = "black",
      fill = "white"
    ) +
    labs(
      x      = "Iteration",
      y      = expression(hat(beta)),
      colour = "Coefficient"
    ) +
    theme_minimal()

  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  p
}
