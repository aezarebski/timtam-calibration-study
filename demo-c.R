targets::tar_source()

uid <- 2
config <- config_c
params <- simulate_parameters(uid, config)
epi <- simulate_epidemic(uid, params, config)
epi_summary <- summarise_epidemic(uid, epi)
epi_msa <- simulate_genomes(epi, params)
mcmc_xml <- get_beast_mcmc_xml(uid, epi_msa, NULL, params, config)
mcmc_samples <- run_beast_mcmc(mcmc_xml, config)
mcmc_samples <- sprintf("out/configuration-c/posterior-samples-%d.log", uid) |>
  read_beast2_log(burn = config$mcmc$num_burn) |>
  dplyr::mutate(log_file = "out/configuration-c/mcmc.xml") |>
  dplyr::select(-Sample) # nolint
mcmc_summary <-
  summarise_post_samples(uid, mcmc_samples, epi_summary, params, config)

png(file.path(dirname(config$files$plot[[1]]), "demo-trace-plot.png"))
mcmc_samples |>
  dplyr::select(posterior, TTR0, TTHistorySizes) |>
  coda::as.mcmc() |>
  plot()
dev.off()

prev_gg <- plot_prevalence_estimates(mcmc_summary)
r0_gg <- plot_r0_estimates(mcmc_summary)
