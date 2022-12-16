#####################
## Configuration A ##
#####################

config_a <- list(
  description = c(
    "Configuration A",
    "Constant rates.",
    "R-naught parameterisation",
    "Assumes unscheduled sequence data and scheduled unsequenced data."
  ),
  cri_coverage = 0.89,
  num_sims = 60,
  param = list(
    r0 = \() runif(n = 1, min = 1.1, max = 1.5),
    sigma = \() 0.2,
    duration = \() 100.0,
    seq_len = \() 1e3,
    seq_rate = \() 1e-3,
    removal_weights = c(4, 1, 3) # for use as an argument to rdirichlet
  ),
  files = list(
    simulation = list(
      remaster_xml = "beast/remaster-a.xml",
      logfile = \(n) sprintf("out/configuration-a/simulation-%d.log", n),
      treefile = \(n) sprintf("out/configuration-a/simulation-%d.tree", n),
      fastafile = \(n) sprintf("out/configuration-a/simulation-%d.fasta", n)
    ),
    mcmc = list(
      template = "beast/timtam-a-template.xml",
      mcmc_xml = \(n) sprintf("out/configuration-a/mcmc-%d.xml", n),
      log_file =
        \(n) sprintf("out/configuration-a/posterior-samples-%d.log", n),
      summary_csv = "out/configuration-a/mcmc-summary.csv"
    ),
    plot = list(
      prevalence_calibration = "out/configuration-a/prevalence-estimates.png",
      r0_calibration = "out/configuration-a/r0-estimates.png"
    )
  ),
  mcmc = list(
    length = 20000000,
    num_burn = 1000,
    min_ess = 200
  ),
  sim_timeout_seconds = 60 * 2,
  mcmc_timeout_seconds = 60 * 60 * 3
)

#####################
## Configuration B ##
#####################

config_b <- config_a
config_b$description <- c(
  "Configuration B",
  "Smaller version of configuration a for testing."
)
config_b$num_sims <- 6
config_b$param$duration <- \() 50
config_b$files <- list(
  simulation = list(
    remaster_xml = "beast/remaster-a.xml",
    logfile = \(n) sprintf("out/configuration-b/simulation-%d.log", n),
    treefile = \(n) sprintf("out/configuration-b/simulation-%d.tree", n),
    fastafile = \(n) sprintf("out/configuration-b/simulation-%d.fasta", n)
  ),
  mcmc = list(
    template = "beast/timtam-a-template.xml",
    mcmc_xml = \(n) sprintf("out/configuration-b/mcmc-%d.xml", n),
    log_file =
      \(n) sprintf("out/configuration-b/posterior-samples-%d.log", n),
    summary_csv = "out/configuration-b/mcmc-summary.csv"
  ),
  plot = list(
    prevalence_calibration = "out/configuration-b/prevalence-estimates.png",
    r0_calibration = "out/configuration-b/r0-estimates.png"
  )
)
config_b$mcmc <- list(
  length = 200000,
  num_burn = 10,
  min_ess = 2
)
config_b$sim_timeout_seconds <- 10
config_b$mcmc_timeout_seconds <- 60

#####################
## Configuration C ##
#####################

config_c <- list(
  description = c(
    "Configuration C",
    "Boom-bust style R0.",
    "R-naught parameterisation",
    "Assumes only unscheduled sequence data"
  ),
  cri_coverage = 0.89,
  num_sims = 6,
  param = list(
    r01 = \() runif(n = 1, min = 1.1, max = 1.5),
    r02 = \() runif(n = 1, min = 0.5, max = 1.0),
    sigma = \() 0.2,
    duration = \() 100.0,
    changeTime = \() 50.0,
    seq_len = \() 1e4,
    seq_rate = \() 1e-2,
    removal_weights = c(4, 1, NA) # for use as an argument to rdirichlet
  ),
  files = list(
    simulation = list(
      remaster_xml = "beast/remaster-c.xml",
      logfile = \(n) sprintf("out/configuration-c/simulation-%d.log", n),
      treefile = \(n) sprintf("out/configuration-c/simulation-%d.tree", n),
      fastafile = \(n) sprintf("out/configuration-c/simulation-%d.fasta", n)
    ),
    mcmc = list(
      template = "beast/timtam-c-template.xml",
      mcmc_xml = \(n) sprintf("out/configuration-c/mcmc-%d.xml", n),
      log_file =
        \(n) sprintf("out/configuration-c/posterior-samples-%d.log", n),
      summary_csv = "out/configuration-c/mcmc-summary.csv"
    ),
    plot = list(
      prevalence_calibration = "out/configuration-c/prevalence-estimates.png",
      r0_calibration = "out/configuration-c/r0-estimates.png"
    )
  ),
  mcmc = list(
    length = 2000000,
    num_burn = 10,
    min_ess = 2
  ),
  sim_timeout_seconds = 10,
  mcmc_timeout_seconds = 5 * 60
)
