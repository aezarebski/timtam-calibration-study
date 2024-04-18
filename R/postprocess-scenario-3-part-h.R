library(ggplot2)

## ================================================================
## The `helper-functions.R` file provides the following functions:
##
## - `read_mcmc`: Read an MCMC log file into a `coda` object.
## - `make_summary`: Summarise an MCMC object with effective size
## - `in_ci`: Check if a value is within the 95% credible interval
## - `ci_width`: Calculate the width of the 95% credible interval
## - `my_ggsave`: Save a ggplot as PNG and SVG in one call
##
source("R/helper-functions.R")
## ================================================================


## The remaster package comes with some helpful tools for working with
## the log-files that it outputs.
source("~/.beast/2.7/remaster/scripts/trajTools.R")

x1 <- loadTrajectories("out/s3/remaster-scenario-3.log")
x2 <- x1 |>
  filter(population == "X") |>
  rename(sample = Sample) |>
  select(t, value, sample)

sim_gg <-
  ggplot(
    data = x2,
    mapping = aes(x = t, y = value)
  ) +
  geom_step(
    mapping = aes(group = sample),
    alpha = 0.6
  ) +
  geom_smooth(
    method = "gam", se = FALSE,
    size = 3,
    colour = "red"
  ) +
  labs(
    x = "Day",
    y = "Prevalence"
  ) +
  theme_bw()

my_ggsave(
  filename = "out/s3/plots/remaster-trajectories-scenario-3.png",
  plot = sim_gg,
  height = 10.5, width = 14.8,
  ## A5 height = 14.8, width = 21.0,
  ## A6 height = 10.5, width = 14.8,
  ## A7 height = 7.4, width = 10.5,
  units = "cm"
)
