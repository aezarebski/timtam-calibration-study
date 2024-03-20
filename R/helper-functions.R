suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coda))


read_mcmc <- function(log_file) {
  if (!file.exists(log_file)) {
    stop("File does not exist: ", log_file)
  }
  log_file |>
    read.table(comment.char = "#", header = TRUE) |>
    select(-Sample) |> # nolint
    as.mcmc()
}
