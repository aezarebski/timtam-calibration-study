config_a <- list(
  description = c(
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


config_b <- config_a
config_b$description <- c("Smaller version of configuration a for testing.")
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
config_b$sim_timeout_seconds = 10
config_b$mcmc_timeout_seconds = 60


#' Sample from Dirichlet distribution.
#'
#' https://en.wikipedia.org/wiki/Dirichlet_distribution#From_gamma_distribution
#'
rdirichlet <- function(shape) {
  tmp <- rgamma(n = 3, shape = c(4, 1, 3), rate = 1)
  tmp / sum(tmp)
}
