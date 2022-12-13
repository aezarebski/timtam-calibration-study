targets::tar_source()

uid <- 1
config <- config_b
params <- simulate_parameters(uid, config)
epi <- simulate_epidemic(uid, params, config)
epi_summary <- summarise_epidemic(uid, epi)
msa <- simulate_genomes(epi, params)
ts_data <- extract_time_series(epi, msa)
mcmc_xml <- get_beast_mcmc_xml(uid, msa, ts_data, params, config)
mcmc_samples <- run_beast_mcmc(mcmc_xml, config)
mcmc_summary <- summarise_post_samples(uid, mcmc_samples, epi_summary, params, config)

png(file.path(dirname(config$files$plot[[1]]), "demo-trace-plot.png"))
mcmc_samples |>
  dplyr::select(posterior, TTR0, TTHistorySizes) |>
  coda::as.mcmc() |>
  plot()
dev.off()
