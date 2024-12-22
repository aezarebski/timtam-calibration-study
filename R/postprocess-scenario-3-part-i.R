library(ggplot2)
library(dplyr)
library(reshape2)


## This data frame contains the total number of sequences and cases
## that were seen in the simulations. I.e. the size of the data set.
total_confirmed_df <-
  "out/s3/total-number-confirmed-cases.csv" |>
  read.csv() |>
  rename(replicate = replicate_num)



.est_prop_summary_df <- function(comb_est_df, var_name,
                                 new_var_name, true_val) {
  comb_est_df |>
    filter(variable == var_name) |>
    mutate(ci_width = fns_5 - fns_1,
           bias = fns_3 - true_val,
           variable = new_var_name) |>
    select(replicate, ci_width, bias, variable)
}

.bind_summaries <- function(summ_df_1, summ_df_2, summ_name) {
  bind_rows(summ_df_1, summ_df_2) |>
    left_join(total_confirmed_df, by = "replicate") |>
    melt(id.vars = c("replicate", "total_confirmed", "variable"),
         variable.name = "property",
         value.name = "value") |>
    mutate(property = factor(property, levels = c("bias", "ci_width")),
           data_source = summ_name)
}

comb_est_3_2_df <-
  "out/s3/summary-combined-estimates-s-3-2.csv" |>
  read.csv()
est_3_2_1_df <- .est_prop_summary_df(comb_est_3_2_df, "R0.1", "Reproduction number 1", 1.850)
est_3_2_2_df <- .est_prop_summary_df(comb_est_3_2_df, "R0.2", "Reproduction number 2", 0.925)
est_3_2_df <- .bind_summaries(est_3_2_1_df, est_3_2_2_df, "Point process")

comb_est_3_3_df <-
  "out/s3/summary-combined-estimates-s-3-3.csv" |>
  read.csv()
est_3_3_1_df <- .est_prop_summary_df(comb_est_3_3_df, "R0.1", "Reproduction number 1", 1.850)
est_3_3_2_df <- .est_prop_summary_df(comb_est_3_3_df, "R0.2", "Reproduction number 2", 0.925)
est_3_3_df <- .bind_summaries(est_3_3_1_df, est_3_3_2_df, "Time series")

plot_df <- bind_rows(est_3_2_df, est_3_3_df)


facet_labels <- c("bias" = "Bias",
                  "ci_width" = "Width of 95% HPD interval")


est_vs_size_gg <-
  ggplot(data = plot_df,
         aes(x = total_confirmed,
             y = value,
             colour = variable)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point() +
  scale_x_log10() +
  scale_colour_manual(values = c("Reproduction number 1" = "#1b9e77",
                                 "Reproduction number 2" = "#7570b3")) +
  labs(x = "Dataset size") +
  facet_grid(property~data_source,
             scales = "free_y",
             labeller = labeller(property = facet_labels)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank())

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
