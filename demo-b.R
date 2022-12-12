targets::tar_source()

uid <- 1
config <- config_b
params <- simulate_parameters(uid, config)
epi <- simulate_epidemic(uid, params, config)
summary <- summarise_epidemic(uid, epi)
msa <- simulate_genomes(uid, epi, params)
