targets::tar_source()

uid <- 1
config <- config_b
params <- simulate_parameters(uid, config)
epi <- simulate_epidemic(uid, params, config)
summary <- summarise_epidemic(uid, epi)
msa <- simulate_genomes(epi, params)
ts_data <- extract_time_series(epi, msa)
mcmc_xml <- get_beast_mcmc_xml(msa, ts_data, params, config)
