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

## ===================================================================
## The `helper-functions.R` file provides the following functions:
##
## - `read_mcmc`: Read an MCMC log file into a `coda` object.
## - `make_summary`: Summarise an MCMC object with effective size
##
source("R/helper-functions.R")
## ===================================================================

summary_gg <- function(plot_df, true_value, name) {
  if (is.null(true_value)) {
    gg <- ggplot() +
      geom_pointrange(
        data = plot_df,
        mapping = aes(
          x = order(sorted_replicate), # nolint
          y = fns_3, # nolint
          ymin = fns_1, # nolint
          ymax = fns_5, # nolint
        )
      ) +
      labs(x = NULL, y = name) +
      theme_bw()
  } else {
    gg <- ggplot() +
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
      theme_bw() +
      theme(legend.position = "top")
  }
  return(gg)
}

## =============================================================================

remaster_node <- "xml/remaster-scenario-3.xml" |> read_xml()
build_node <- read_xml("build.xml")
num_replicates <- build_node |>
  xml_find_first(xpath = "//property[@name='numSims']") |>
  xml_attr("value") |>
  as.integer()

summary_list <-
  map(
    seq.int(num_replicates),
    \(ix) make_summary(
            read_mcmc(
              str_interp("out/s3/timtam-scenario-3-3-sample-$[03d]{ix}.log")
            ),
            ix)
  )

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
            file = "out/s3/plots/effective-sample-sizes-s-3-3.csv",
            sep = ",",
            row.names = FALSE)

ggsave(filename = "out/s3/plots/effective-sample-sizes-s-3-3.png",
       plot = ess_gg,
       height = 10.5, width = 14.8,
       units = "cm")

## =============================================================================

birth_rate_df <- comb_est_df |>
  filter(variable == "TTBirthRate.1") |>
  mutate(sorted_replicate = order(fns_3))

birth_rate_gg <- summary_gg(birth_rate_df, 0.185, "Birth rate")

ggsave(
  filename = "out/s3/plots/birth-rate-1-estimates-s-3-3.png",
  plot = birth_rate_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

birth_rate_df <- comb_est_df |>
  filter(variable == "TTBirthRate.2") |>
  mutate(sorted_replicate = order(fns_3))

birth_rate_gg <- summary_gg(birth_rate_df, 0.0925, "Birth rate")

ggsave(
  filename = "out/s3/plots/birth-rate-2-estimates-s-3-3.png",
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
  filename = "out/s3/plots/sampling-rate-estimates-s-3-3.png",
  plot = sampling_rate_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

## =============================================================================

nu_prob_df <- comb_est_df |>
  filter(variable == "TTNuProb") |>
  mutate(sorted_replicate = order(fns_3))

nu_prob_gg <- summary_gg(nu_prob_df, NULL, "Nu probability")

ggsave(
  filename = "out/s3/plots/nu-prob-estimates-s-3-3.png",
  plot = nu_prob_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

## =============================================================================
