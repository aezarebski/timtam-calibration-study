# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tibble"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

list(
  tar_target(uid_b, seq.int(config_b$num_sims)),
  tar_target(
    params_b,
    simulate_parameters(uid_b, config_b),
    pattern = map(uid_b)
  ),
  tar_target(
    epi_b,
    simulate_epidemic(uid_b, params_b, config_b),
    pattern = map(uid_b, params_b)
  ),
  tar_target(
    epi_summary_b,
    summarise_epidemic(uid_b, epi_b),
    pattern = map(uid_b, epi_b)
  ),
  tar_target(
    epi_msa_b,
    simulate_genomes(epi_b, params_b),
    pattern = map(epi_b, params_b)
  ),
  tar_target(
    epi_ts_data_b,
    extract_time_series(epi_b, epi_msa_b),
    pattern = map(epi_b, epi_msa_b)
  ),
  tar_target(
    mcmc_xml_b,
    get_beast_mcmc_xml(uid_b, epi_msa_b, epi_ts_data_b, params_b, config_b),
    pattern = map(uid_b, epi_msa_b, epi_ts_data_b, params_b)
  ),
  tar_target(
    mcmc_samples_b,
    run_beast_mcmc(mcmc_xml_b, config_b),
    pattern = map(mcmc_xml_b)
  ),
  tar_target(
    mcmc_summary_b,
    summarise_post_samples(
      uid_b, mcmc_samples_b, epi_summary_b, params_b, config_b
    ),
    pattern = map(uid_b, mcmc_samples_b, epi_summary_b, params_b)
  ),
  tar_target(
    prev_ests_gg_b,
    plot_prevalence_estimates(mcmc_summary_b)
  ),
  tar_target(
    r0_ests_gg_b,
    plot_r0_estimates(mcmc_summary_b)
  ),
  tar_target(
    comb_gg_b,
    combined_figure(prev_ests_gg_b, r0_ests_gg_b))
)
