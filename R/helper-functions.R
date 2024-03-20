suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coda))


#' Read MCMC log file
#'
#' Read MCMC log file and return an MCMC object
#'
#' @param log_file Path to the log file
#' @return MCMC object
#'
read_mcmc <- function(log_file) {
  if (!file.exists(log_file)) {
    stop("File does not exist: ", log_file)
  }
  log_file |>
    read.table(comment.char = "#", header = TRUE) |>
    select(-Sample) |> # nolint
    as.mcmc()
}
