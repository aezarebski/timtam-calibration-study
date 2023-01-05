suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(reshape2))
set.seed(1)

ix <- 1

read_mcmc <- function(ix) {
  str_interp("out/s2/timtam-scenario-2-sample-$[03d]{ix}.log") |>
    read.table(comment.char = "#", header = T) |>
    select(-Sample) |> # nolint
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

summary_gg <- function(plot_df, true_value, name) {
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
    geom_hline(yintercept = true_value, linetype = "dashed") +
    labs(x = NULL, y = name, shape = "Contains true value") +
    scale_shape_manual(breaks = c(FALSE, TRUE), values = c(1, 16)) +
    scale_linetype_manual(breaks = c(FALSE, TRUE), values = c("dashed", "solid")) +
    theme_bw() +
    theme(legend.position = "top")
}

## =============================================================================

remaster_node <- "remaster-scenario-2.xml" |> read_xml()
num_replicates <- remaster_node |>
  xml_find_first("//run") |>
  xml_attr("nSims") |>
  as.integer()

summary_list <- map(seq.int(num_replicates), \(x) make_summary(read_mcmc(x), x))

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

ggsave(filename = "out/s2/plots/effective-sample-sizes.png",
       plot = ess_gg,
       height = 10.5, width = 14.8,
       units = "cm")

## =============================================================================

birth_rate_df <- comb_est_df |>
  filter(variable == "TTBirthRate") |>
  mutate(sorted_replicate = order(fns_3))

birth_rate_gg <- summary_gg(birth_rate_df, 0.185, "Birth rate")

ggsave(
  filename = "out/s2/plots/birth-rate-estimates.png",
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
  filename = "out/s2/plots/sampling-rate-estimates.png",
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
  filename = "out/s2/plots/omega-rate-estimates.png",
  plot = omega_rate_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

## =============================================================================
