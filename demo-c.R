## targets::tar_source()
source("R/configuration.R")
source("R/functions.R")
uid <- 2
config <- config_c
params <- simulate_parameters(uid, config)
epi <- simulate_epidemic(uid, params, config)
epi_summary <- summarise_epidemic(uid, epi)
print(epi_summary)
epi_msa <- simulate_genomes(epi, params)
## Save a copy of the sequences so that they can be used in a BDSKY analysis.
phangorn::write.phyDat(
            epi_msa,
            config$files$simulation$fastafile(uid),
            format = "fasta"
          )
mcmc_xml <- get_beast_mcmc_xml(uid, epi_msa, NULL, params, config)
xml2::write_xml(mcmc_xml, "out/configuration-c/demo-c-mcmc.xml")
mcmc_samples <- run_beast_mcmc(mcmc_xml, config)
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
