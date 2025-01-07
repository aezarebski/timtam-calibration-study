## This script generates a plot that shows how the error and width of
## the 95% HPD interval of the estimates of the reproduction numbers
## change as the size of the data set increases.
##
## Technical note
## --------------
##
## We use the "error" rather than the "bias" of the estimates. We
## compute the error as the value of the point estimate minus the true
## value.
##
library(ggplot2)
library(dplyr)
library(reshape2)
DRY_RUN <- FALSE # Print the plot if TRUE, save it if FALSE.

## True parameter values used in the simulation. Yes, it's gross that
## they are hard-coded, but at least they are defined in one place
## rather than appearing as magic numbers scattered throughout the
## code.
true_r0_1 <- 1.850
true_r0_2 <- 0.925

## This data frame contains the total number of sequences and cases
## that were seen in the simulations. I.e. the size of the data set.
## Calculating this from scratch is painful so better to just read in
## the saved version that already exists.
total_confirmed_df <-
  "out/s3/total-number-confirmed-cases.csv" |>
  read.csv() |>
  rename(replicate = replicate_num)


.est_prop_summary_df <- function(comb_est_df, var_name,
                                 new_var_name, true_val) {
  comb_est_df |>
    filter(variable == var_name) |>
    mutate(ci_width = fns_5 - fns_1,
           error = fns_3 - true_val,
           variable = new_var_name) |>
    select(replicate, ci_width, error, variable)
}

.bind_summaries <- function(summ_df_1, summ_df_2, summ_name) {
  bind_rows(summ_df_1, summ_df_2) |>
    left_join(total_confirmed_df, by = "replicate") |>
    melt(id.vars = c("replicate", "total_confirmed", "variable"),
         variable.name = "property",
         value.name = "value") |>
    mutate(property = factor(property, levels = c("error", "ci_width")),
           data_source = summ_name)
}

comb_est_3_2_df <-
  "out/s3/summary-combined-estimates-s-3-2.csv" |>
  read.csv()
est_3_2_1_df <- .est_prop_summary_df(comb_est_3_2_df, "R0.1", "Reproduction number 1", true_r0_1)
est_3_2_2_df <- .est_prop_summary_df(comb_est_3_2_df, "R0.2", "Reproduction number 2", true_r0_2)
est_3_2_df <- .bind_summaries(est_3_2_1_df, est_3_2_2_df, "Point process")

comb_est_3_3_df <-
  "out/s3/summary-combined-estimates-s-3-3.csv" |>
  read.csv()
est_3_3_1_df <- .est_prop_summary_df(comb_est_3_3_df, "R0.1", "Reproduction number 1", true_r0_1)
est_3_3_2_df <- .est_prop_summary_df(comb_est_3_3_df, "R0.2", "Reproduction number 2", true_r0_2)
est_3_3_df <- .bind_summaries(est_3_3_1_df, est_3_3_2_df, "Time series")

plot_df <- bind_rows(est_3_2_df, est_3_3_df)


facet_labels <- c("error" = "Error of the point estimate",
                  "ci_width" = "Width of 95% HPD interval")


## The `geom_smooth` layer is used to show average error across the
## set of replicates so actually measures the bias. But since each
## individual point is the error, this is what we have labelled it as.
est_vs_size_gg <-
  ggplot(data = plot_df,
         aes(x = total_confirmed,
             y = value,
             colour = variable)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point(size = 1.5) +
  geom_smooth(data = subset(plot_df, property == "error"),
              method = "lm", formula = y ~ 1, se = TRUE,
              show.legend = FALSE) +
  scale_x_log10() +
  scale_colour_manual(values = c("Reproduction number 1" = "#1b9e77",
                                 "Reproduction number 2" = "#7570b3"),
                            name = "Parameter") +
  labs(x = "Dataset size") +
  facet_grid(property~data_source,
             scales = "free_y",
             labeller = labeller(property = facet_labels)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "top")

if (DRY_RUN) {
  print(est_vs_size_gg)
} else {
  fig_width <- 6
  fig_height <- 6
  ggsave("out/s3/plots/summary-est-vs-size.png",
         plot = est_vs_size_gg,
         width = fig_width,
         height = fig_height,
         dpi = 300)
  ggsave("out/s3/plots/summary-est-vs-size.svg",
         plot = est_vs_size_gg,
         width = fig_width,
         height = fig_height)
}
