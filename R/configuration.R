config_a <- list(
  description = "Constant rates.",
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
      logfile = \(n) sprintf("out/configuration-a/simulation-%d.log", n),
      treefile = \(n) sprintf("out/configuration-a/simulation-%d.tree", n),
      fastafile = \(n) sprintf("out/configuration-a/simulation-%d.fasta", n)
    ),
    mcmc = list(
      template = "beast/timtam-a-template.xml",
      mcmc_xml = \(n) sprintf("out/configuration-a/mcmc-%d.xml", n),
      summary_csv = "out/mcmc-summary.csv"
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


#' Sample from Dirichlet distribution.
#'
#' https://en.wikipedia.org/wiki/Dirichlet_distribution#From_gamma_distribution
#'
rdirichlet <- function(shape) {
  tmp <- rgamma(n = 3, shape = c(4, 1, 3), rate = 1)
  tmp / sum(tmp)
}
