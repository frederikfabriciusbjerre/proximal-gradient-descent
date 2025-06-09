library(tidyverse)

# ──────────────────────────────── Plot theme ──────────────────────────────────
small_plot_theme <- function() {
  theme_minimal() %+replace% # replace elements we want to change
    theme(
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      plot.title = element_text(size = 8)
    )
}
# ──────────────────────────────────────────────────────────────────────────────

# read csv
results_df <- read.csv("results/pcd_sim_results_1_10.csv")
results_df_mcp <- read.csv("results/pcd_sim_results_mcp_1_10.csv")
metric_cols <- c("l1", "mse_full", "intercept_err..Intercept.")

# function that throws out scenarios that not all models
keep_complete_scenarios <- function(df,
                                    strategy_col = "intercept_update") {
  drop_cols <- c(
    strategy_col,
    grep("(_iter$|_err|mse|l1)", names(df), value = TRUE)
  )
  scenario_cols <- setdiff(names(df), drop_cols)

  n_strats <- df |>
    distinct(.data[[strategy_col]]) |>
    nrow()

  df |>
    group_by(across(all_of(scenario_cols))) |>
    filter(n_distinct(.data[[strategy_col]]) == n_strats) |>
    ungroup()
}

conv_summary <- function(df, strategy_col = "intercept_update") {
  df |>
    group_by(.data[[strategy_col]]) |>
    summarise(
      n_runs = n(),
      mean = mean(conv_iter, na.rm = TRUE),
      sd = sd(conv_iter, na.rm = TRUE),
      median = median(conv_iter, na.rm = TRUE),
      IQR = IQR(conv_iter, na.rm = TRUE),
      .groups = "drop"
    )
}

# retain only complete scenarios
results_df_binomial <- keep_complete_scenarios(results_df)
results_df_mcp_binomial <- keep_complete_scenarios(results_df_mcp)

summary_results <- conv_summary(results_df)
summary_results_mcp <- conv_summary(results_df_mcp)

# concatenate results
results_df_combined <- bind_rows(
  results_df_binomial |> mutate(prox_fun = "soft"),
  results_df |> filter(family == "gaussian") |> mutate(prox_fun = "soft"),
  results_df_mcp_binomial |> mutate(prox_fun = "mcp"),
  results_df_mcp |> filter(family == "gaussian") |> mutate(prox_fun = "mcp")
)

# ──────────────────────────── Convergence plots ───────────────────────────────
p <- ggplot(
  results_df_binomial,
  aes(x = intercept_update, y = conv_iter, fill = method)
) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Sweeps until convergence for binomial family",
    x = "Intercept-update strategy",
    y = "converged sweeps (log-scale)"
  )

ggsave(p,
  filename = "plots/soft_binom_convergence.pdf",
  width = 6, height = 4
)

p <- ggplot(
  results_df |> filter(family == "gaussian"),
  aes(x = intercept_update, y = conv_iter, fill = method)
) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Sweeps until convergence for Gaussian family",
    x = "Intercept-update strategy",
    y = "converged sweeps (log-scale)"
  )

ggsave(p,
  filename = "plots/soft_gauss_convergence.pdf",
  width = 6, height = 4
)

# facet grid over family
p <- ggplot(
  results_df,
  aes(x = intercept_update, y = conv_iter, fill = method)
) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(family ~ ., scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Sweeps until convergence",
    x = "Intercept-update strategy",
    y = "converged sweeps (log-scale)"
  )

ggsave(p,
  filename = "plots/soft_both_convergence.pdf",
  width = 6, height = 8
)

# ───────────────────────────────── L1 plots ───────────────────────────────────
p <- ggplot(
  results_df,
  aes(x = intercept_update, y = l1, fill = method)
) +
  geom_boxplot() +
  facet_grid(family ~ ., scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "L1 distance of ß_hat vs ß",
    x = "Intercept-update strategy",
    y = "L1 error"
  )
ggsave(p,
  filename = "plots/soft_both_l1.pdf",
  width = 9, height = 8
)
# stratified over design for binomial
p <- ggplot(
  results_df_binomial,
  aes(x = intercept_update, y = l1, fill = method)
) +
  geom_boxplot() +
  facet_grid(. ~ design) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "L1 distance of ß_hat vs ß",
    x = "Design",
    y = "L1 error"
  )
ggsave(p,
  filename = "plots/soft_binomial_l1_design_stratified.pdf",
  width = 9, height = 4
)
# stratified over design for gaussian
p <- ggplot(
  results_df |> filter(family == "gaussian"),
  aes(x = intercept_update, y = l1, fill = method)
) +
  geom_boxplot() +
  facet_grid(. ~ design) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "L1 distance of ß_hat vs ß",
    x = "Design",
    y = "L1 error"
  )
ggsave(p,
  filename = "plots/soft_gauss_l1_design_stratified.pdf",
  width = 9, height = 4
)
# ────────────────────────────── Intercept plots ───────────────────────────────
p <- ggplot(
  results_df,
  aes(x = intercept_update, y = intercept_err..Intercept., fill = method)
) +
  geom_boxplot() +
  facet_grid(family ~ ., scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Intercept error",
    x = "Intercept-update strategy",
    y = "Intercept error"
  )
ggsave(p,
  filename = "plots/soft_both_intercept.pdf",
  width = 9, height = 8
)

# stratified over design for binomial
p <- ggplot(
  results_df_binomial,
  aes(x = intercept_update, y = intercept_err..Intercept., fill = method)
) +
  geom_boxplot() +
  facet_grid(prox_fun ~ design) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Intercept error for binomial family",
    x = "Intercept-update strategy",
    y = "Intercept error"
  )
ggsave(p,
  filename = "plots/soft_binomial_intercept_l1_design_stratified.pdf",
  width = 9, height = 4
)

# stratified over design for gaussian
p <- ggplot(
  results_df |> filter(family == "gaussian"),
  aes(x = intercept_update, y = intercept_err..Intercept., fill = method)
) +
  geom_boxplot() +
  coord_cartesian(ylim = quantile(results_df$intercept_err..Intercept., c(0, 1))) +
  facet_grid(prox_fun ~ design, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Intercept error for Gaussian family",
    x = "Intercept-update strategy",
    y = "Intercept error"
  )
ggsave(p,
  filename = "plots/soft_gauss_intercept_l1_design_stratified.pdf",
  width = 9, height = 4
)

# ──────────────────────────────────── MCP ─────────────────────────────────────
p <- ggplot(
  {
    df <- results_df_combined

    # threshold that marks the highest 0.1 % in the mcp–binomial subset
    thr <- quantile(
      df$l1[df$prox_fun == "mcp" & df$family == "binomial"],
      0.999
    ) # keep 99.9 %

    df %>%
      filter(!(prox_fun == "mcp" & family == "binomial" & l1 > thr)) # drop outliers
  } |> filter(family == "binomial"),
  aes(x = intercept_update, y = l1, fill = method)
) +
  geom_boxplot() +
  facet_grid(prox_fun ~ family, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "L1 distance of ß_hat vs ß for MCP vs Soft",
    x = "Intercept-update strategy",
    y = "L1 error"
  )

ggsave(p,
  filename = "plots/mcp_binomial_l1.pdf",
  width = 9, height = 8
)

# convergence comparison
ggplot(
  results_df_combined,
  aes(x = intercept_update, y = conv_iter, fill = method)
) +
  geom_boxplot() +
  facet_grid(prox_fun ~ family, scales = "free_y") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Sweeps until convergence for binomial family (MCP)",
    x = "Intercept-update strategy",
    y = "converged sweeps (log-scale)"
  )
