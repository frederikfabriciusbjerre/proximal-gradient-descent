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

# ──────────────────────────────── Read data ───────────────────────────────────

results_df <- read.csv("results/full_sim_results.csv",
  header = TRUE,
  row.names = 1
)
results_df_mcp <- filter(results_df, prox_fun == "mcp")
results_df_soft <- filter(results_df, prox_fun == "soft")
metric_cols <- c("l1", "mse_full", "intercept_err..Intercept.")


# ──────────────────────────────── Clean data ──────────────────────────────────

# function that throws out scenarios that not all models
keep_complete_scenarios <- function(df) {
  strategy_col <- "intercept_update"

  drop_cols <- c(
    strategy_col,
    "n_updates",
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

# retain only complete scenarios
results_df_soft <- keep_complete_scenarios(results_df_soft)
results_df_mcp <- keep_complete_scenarios(results_df_mcp)

# ──────────────────────────────── Summaries ───────────────────────────────────

conv_summary <- function(df, strategy_col = "intercept_update") {
  df |>
    group_by(.data[[strategy_col]]) |>
    summarise(
      "Mean" = mean(conv_iter, na.rm = TRUE),
      "sd" = sd(conv_iter, na.rm = TRUE),
      "Median" = median(conv_iter, na.rm = TRUE),
      "IQR" = IQR(conv_iter, na.rm = TRUE),
      .groups = "drop"
    )
}
n_update_summary <- function(df, strategy_col = "intercept_update") {
  df |>
    group_by(.data[[strategy_col]]) |>
    summarise(
      "Mean" = mean(n_updates, na.rm = TRUE),
      "sd" = sd(n_updates, na.rm = TRUE),
      "Median" = median(n_updates, na.rm = TRUE),
      "IQR" = IQR(n_updates, na.rm = TRUE),
      .groups = "drop"
    )
}

conv_results_soft <- conv_summary(results_df_soft)
conv_results_mcp <- conv_summary(results_df_mcp)
print(conv_results_soft)
print(conv_results_mcp)

n_update_results_soft <- n_update_summary(results_df_soft)
n_update_results_mcp <- n_update_summary(results_df_mcp)
print(n_update_results_soft)
print(n_update_results_mcp)

# ──────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────── Plots ────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────── Convergence plots ───────────────────────────────

p <- ggplot(
  results_df_soft |> filter(family == "binomial"),
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
  results_df_soft |> filter(family == "gaussian") |> filter(method == "gradient"),
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
  results_df_soft,
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

# ───────────────────────────── n updates plots ────────────────────────────────

p <- ggplot(
  results_df_soft |> filter(family == "binomial"),
  aes(x = intercept_update, y = n_updates, fill = method)
) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Number of updates until convergence for binomial family",
    x = "Intercept-update strategy",
    y = "converged sweeps (log-scale)"
  )

ggsave(p,
  filename = "plots/soft_binom_n_updates.pdf",
  width = 6, height = 4
)

p <- ggplot(
  results_df_soft |> filter(family == "gaussian") |> filter(method == "gradient"),
  aes(x = intercept_update, y = n_updates, fill = method)
) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Number of updates until convergence for binomial family",
    x = "Intercept-update strategy",
    y = "converged sweeps (log-scale)"
  )

ggsave(p,
  filename = "plots/soft_gaussian_n_updates.pdf",
  width = 6, height = 4
)

# ───────────────────────────────── L1 plots ───────────────────────────────────

p <- ggplot(
  results_df_soft,
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
  results_df_soft |> filter(family == "binomial"),
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
  results_df_soft |> filter(family == "gaussian") |>
    filter(method == "gradient"),
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
  results_df_soft,
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
  results_df_soft |> filter(family == "binomial"),
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
  results_df_soft |> filter(family == "gaussian") |>
    filter(method == "gradient"),
  aes(x = intercept_update, y = intercept_err..Intercept., fill = method)
) +
  geom_boxplot() +
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
  results_df |> filter(family == "binomial"),
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
  results_df,
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
