library(ggplot2)

## The remaster package comes with some helpful tools for working with
## the log-files that it outputs.
source("~/.beast/2.7/remaster/scripts/trajTools.R")

x1 <- loadTrajectories("out/s3/remaster-scenario-3.log")
x2 <- x1 |>
  filter(population == "X") |>
  rename(sample = Sample) |>
  select(t, value, sample)

sim_gg <-
  ggplot(data = x2,
         mapping = aes(x = t, y = value)) +
  geom_step(mapping = aes(group = sample),
            alpha = 0.6) +
  geom_smooth(method = "gam", se = FALSE,
              size = 3,
              colour = "red") +
  labs(
    x = "Day",
    y = "Prevalence"
  ) +
  theme_bw()

ggsave(filename = "out/s3/plots/remaster-trajectories-scenario-3.png",
       plot = sim_gg,
       height = 10.5, width = 14.8,
       ## A5 height = 14.8, width = 21.0,
       ## A6 height = 10.5, width = 14.8,
       ## A7 height = 7.4, width = 10.5,
       units = "cm")
