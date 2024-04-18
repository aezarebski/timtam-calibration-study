library(ggplot2)
library(dplyr)

final_prevalence_df <-
  "out/s3/final-simulation-state.csv" |>
  read.csv() |>
  select(X, sample) |>
  mutate(prevalence = X, replicate_num = sample + 1)

total_confirmed_df <-
  "out/s3/total-number-confirmed-cases.csv" |>
  read.csv()

plot_df <- total_confirmed_df |>
  left_join(final_prevalence_df, by = "replicate_num") |>
  select(total_confirmed, prevalence)

prev_vs_num_gg <- ggplot() +
  geom_point(data = plot_df, mapping = aes(x = total_confirmed, y = prevalence)) +
  labs(
    x = "Total number of confirmed cases",
    y = "Final prevalence of infection"
  ) +
  theme_bw()

ggsave(
  filename = "out/s3/plots/prevalence-data-set-size-plot.png",
  plot = prev_vs_num_gg,
  height = 10.5, width = 10.5,
  units = "cm"
)
