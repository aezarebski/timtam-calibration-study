suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(reshape2))
set.seed(1)

occurrence_files <- list.files("out/s3", pattern = ".*ssv", full.names = TRUE)
sequence_files <- list.files("out/s3", pattern = ".*fasta", full.names = TRUE)

replicate_number <- function(fn) {
  fn |>
    basename() |>
    str_match("([0-9]+)\\.[a-z]+$") |>
    pluck(2) |>
    as.integer()
}

read_num_sequences <- function(seq_file) {
  num_seq <- seq_file |>
    read.FASTA() |>
    length()
  data.frame(num_seq = num_seq, replicate_num = replicate_number(seq_file))
}

read_num_occurrences <- function(occ_file) {
  num_occ <- occ_file |>
    readLines() |>
    str_split(" ") |>
    pluck(1) |>
    length()
  data.frame(num_occ = num_occ, replicate_num = replicate_number(occ_file))
}

seq_df <- sequence_files |>
  map(read_num_sequences) |>
  bind_rows()
occ_df <- occurrence_files |>
  map(read_num_occurrences) |>
  bind_rows()

data_size_df <- seq_df |>
  left_join(occ_df, by = "replicate_num") |>
  melt(id.vars = "replicate_num")

tmp_1 <- data_size_df |>
  group_by(variable) |>
  summarise(mean = mean(value)) |>
  as.data.frame()

data_size_1_gg <-
  ggplot() +
  geom_point(
    data = data_size_df,
    mapping = aes(
      x = replicate_num,
      y = value,
      colour = variable
    )
  ) +
  geom_hline(
    data = tmp_1,
    mapping = aes(
      yintercept = mean,
      colour = variable
    ),
    linetype = "dashed"
  ) +
  theme_bw()

ggsave(
  filename = "out/s3/plots/dataset-size-1.png",
  plot = data_size_1_gg,
  height = 10.5, width = 14.8,
  units = "cm"
)

data_size_2_gg <-
  ggplot() +
  geom_histogram(
    data = data_size_df,
    mapping = aes(
      x = value,
      y = after_stat(density),
      fill = variable,
      colour = variable
    ),
    alpha = 0.2,
    bins = 30
  ) +
  geom_vline(
    data = tmp_1,
    mapping = aes(
      xintercept = mean,
      colour = variable
    ),
    linetype = "dashed"
  ) +
  facet_wrap(~variable) +
  theme_bw()

ggsave(
  filename = "out/s3/plots/dataset-size-2.png",
  plot = data_size_2_gg,
  height = 14.8, width = 21.0,
  units = "cm"
)

total_confirmed_df <- data_size_df |>
  group_by(replicate_num) |>
  summarise(total_confirmed = sum(value)) |>
  as.data.frame()

write.table(
  x = total_confirmed_df,
  file = "out/s3/total-number-confirmed-cases.csv",
  sep = ",",
  row.names = FALSE
)
