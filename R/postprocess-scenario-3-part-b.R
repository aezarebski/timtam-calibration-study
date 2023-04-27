suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(reshape2))
set.seed(1)

## library(cowplot)
## library(RColorBrewer)
## palette <- brewer.pal(5, "Dark2")[1:3]
highlight_col_hex <- "#1B9E77"
ix <- 1

remaster_xml <- "xml/remaster-scenario-3.xml"
if (!file.exists(remaster_xml)) {
  stop("File does not exist: ", remaster_xml)
}

read_mcmc <- function(log_file) {
  ## this should throw a helpful error if the file doesn't exist
  if (!file.exists(log_file)) {
    stop("File does not exist: ", log_file)
  }
  log_file |>
    read.table(comment.char = "#", header = T) |>
    select(-Sample) |> # nolint
    mutate(
      R0.1 = TTBirthRate.1 / (TTSamplingRate + TTOmegaRate + 0.046),
      R0.2 = TTBirthRate.2 / (TTSamplingRate + TTOmegaRate + 0.046)
    ) |>
    as.mcmc()
}

make_summary <- function(mcmc_obj, ix) {
  summary_df <- summary(mcmc_obj)$quantiles
  tmp <- rownames(summary_df)
  summary_df <- summary_df |>
    as.data.frame() |>
    mutate(variable = tmp, replicate = ix) |>
    rename(
      fns_1 = "2.5%",
      fns_2 = "25%",
      fns_3 = "50%",
      fns_4 = "75%",
      fns_5 = "97.5%",
    )

  eff_size_df <- mcmc_obj |>
    effectiveSize() |>
    t() |>
    as.data.frame() |>
    mutate(replicate = ix)

  return(list(
    summary = summary_df,
    effective_sizes = eff_size_df
  ))
}

summary_gg <- function(plot_df, true_value, name, hline_col="black") {
  ggplot() +
    geom_pointrange(
      data = plot_df,
      mapping = aes(
        x = order(sorted_replicate), # nolint
        y = fns_3, # nolint
        ymin = fns_1, # nolint
        ymax = fns_5, # nolint
        shape = fns_1 < true_value & true_value < fns_5
      )
    ) +
    geom_hline(yintercept = true_value, linetype = "dashed", colour = hline_col, linewidth = 1.0) +
    labs(x = NULL, y = name, shape = "Contains true value") +
    scale_shape_manual(breaks = c(FALSE, TRUE), values = c(1, 16)) +
    scale_linetype_manual(breaks = c(FALSE, TRUE), values = c("dashed", "solid")) +
    theme_bw() +
    theme(legend.position = "top")
}

## =============================================================================

remaster_node <- remaster_xml |> read_xml()
build_node <- read_xml("build.xml")
num_replicates <- build_node |>
  xml_find_first(xpath = "//property[@name='numSims']") |>
  xml_attr("value") |>
  as.integer()

summary_list <- map(seq.int(num_replicates), \(ix) make_summary(read_mcmc(str_interp("out/s3/timtam-scenario-3-2-sample-$[03d]{ix}.log")), ix))

comb_est_df <- summary_list |>
  map(\(x) x$summary) |>
  bind_rows()

comb_effsize_df <- summary_list |>
  map(\(x) x$effective_sizes) |>
  bind_rows() |>
  melt(id.vars = "replicate")

## =============================================================================

ess_gg <- ggplot() +
  geom_point(
    data = comb_effsize_df,
    mapping = aes(x = replicate, y = value, colour = variable)
  ) +
  geom_hline(yintercept = 200, linetype = "dashed") +
  labs(x = NULL, y = "Effective sample size") +
  theme_bw()

write.table(x = comb_effsize_df,
            file = "out/s3/plots/effective-sample-sizes-s-3-2.csv",
            sep = ",",
            row.names = FALSE)

ggsave(filename = "out/s3/plots/effective-sample-sizes-s-3-2.png",
       plot = ess_gg,
       height = 10.5, width = 14.8,
       units = "cm")

## =============================================================================

birth_rate_df <- comb_est_df |>
  filter(variable == "TTBirthRate.1") |>
  mutate(sorted_replicate = order(fns_3))

birth_rate_gg <- summary_gg(birth_rate_df, 0.185, "Birth rate")

ggsave(
  filename = "out/s3/plots/birth-rate-1-estimates-s-3-2.png",
  plot = birth_rate_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

birth_rate_df <- comb_est_df |>
  filter(variable == "TTBirthRate.2") |>
  mutate(sorted_replicate = order(fns_3))

birth_rate_gg <- summary_gg(birth_rate_df, 0.0925, "Birth rate")

ggsave(
  filename = "out/s3/plots/birth-rate-2-estimates-s-3-2.png",
  plot = birth_rate_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

## =============================================================================

sampling_rate_df <- comb_est_df |>
  filter(variable == "TTSamplingRate") |>
  mutate(sorted_replicate = order(fns_3))

sampling_rate_gg <- summary_gg(sampling_rate_df, 0.008, "Sampling rate")

ggsave(
  filename = "out/s3/plots/sampling-rate-estimates-s-3-2.png",
  plot = sampling_rate_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

## =============================================================================

omega_rate_df <- comb_est_df |>
  filter(variable == "TTOmegaRate") |>
  mutate(sorted_replicate = order(fns_3))

omega_rate_gg <- summary_gg(omega_rate_df, 0.046, "Omega rate")

ggsave(
  filename = "out/s3/plots/omega-rate-estimates-s-3-2.png",
  plot = omega_rate_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

## =============================================================================

final_prevalence_df <-
  "out/s3/final-simulation-state.csv" |>
  read.csv() |>
  select(X, sample) |>
  mutate(prevalence = X, replicate = sample + 1)

history_sizes_df <- comb_est_df |>
  filter(variable == "HistorySizes") |>
  inner_join(final_prevalence_df, by = "replicate") |>
  mutate(sorted_replicate = order(prevalence)) |>
  mutate(contains_truth = fns_1 < prevalence & prevalence < fns_5)

dummy_data <- data.frame(
  x = NaN, y = NaN,
  contains_truth = c(FALSE, TRUE)
)

history_size_gg <-
  ggplot(data = history_sizes_df) +
  geom_pointrange(
    mapping = aes(
      x = order(sorted_replicate), # nolint
      y = fns_3, # nolint
      ymin = fns_1, # nolint
      ymax = fns_5, # nolint
      shape = contains_truth,
      colour = "estimate"
    )
  ) +
  geom_point(
    mapping = aes(
      x = order(sorted_replicate),
      y = prevalence,
      colour = "truth"
    ),
    size = 2.5
  ) +
  geom_point(
    data = dummy_data,
    mapping = aes(x = x, y = y, shape = contains_truth),
    inherit.aes = FALSE,
    na.rm = TRUE
  ) +
  scale_y_log10() +
  labs(x = NULL, y = "Prevalence", shape = "Contains true value") +
  scale_colour_manual("Type",
    breaks = c("estimate", "truth"), values = c("black", highlight_col_hex)) +
  scale_shape_manual(breaks = c(FALSE, TRUE),
                     values = c(1, 16),
                     labels = c("No", "Yes"),
                     drop = FALSE) +
  scale_linetype_manual(breaks = c(FALSE, TRUE), values = c("dashed", "solid")) +
  theme_bw() +
  theme(legend.position = "top")

## =============================================================================

r0_1_df <- comb_est_df |>
  filter(variable == "R0.1")
r0_1_df$sorted_replicate <- history_sizes_df$sorted_replicate

r0_1_gg <- summary_gg(r0_1_df, 1.85, "Reproduction number 1", hline_col=highlight_col_hex) +
  theme(legend.position = "NONE")

## =============================================================================

r0_2_df <- comb_est_df |>
  filter(variable == "R0.2") |>
  mutate(sorted_replicate = order(fns_3))
r0_2_df$sorted_replicate <- history_sizes_df$sorted_replicate

r0_2_gg <- summary_gg(r0_2_df, 0.925, "Reproduction number 2", hline_col=highlight_col_hex) +
  theme(legend.position = "NONE")


## =============================================================================

combined_gg <- plot_grid(
  history_size_gg + theme(axis.text.x = element_blank()),
  r0_1_gg + theme(axis.text.x = element_blank()),
  r0_2_gg,
  ncol = 1,
  align = "v"
)

ggsave(
  filename = "out/s3/plots/combined-r0-prevalence-estimates-s-3-2.png",
  plot = combined_gg,
  height = 3 * 10.5, width = 14.8,
  units = "cm"
)
