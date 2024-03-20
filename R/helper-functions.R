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


#' Make summary of MCMC object
#'
#' Make summary of MCMC object and return a list of data frames
#' with summary statistics and effective sizes
#'
#' @param mcmc_obj MCMC object
#' @param ix Replicate index
#' @return List of data frames
#'
make_summary <- function(mcmc_obj, ix) {
  if (!is(mcmc_obj, "mcmc")) {
    stop("Input is not an MCMC object")
  }
  summary_df <- summary(mcmc_obj)$quantiles
  tmp <- rownames(summary_df)
  summary_df <- summary_df |>
    as.data.frame() |>
    mutate(variable = tmp, replicate = ix) |>
    rename(
      fns_1 = "2.5%",
      fns_2 = "25%",
      fns_3 = "50%",
      fns_4 = "75%",
      fns_5 = "97.5%",
    )

  eff_size_df <- mcmc_obj |>
    effectiveSize() |>
    t() |>
    as.data.frame() |>
    mutate(replicate = ix)

  list(
    summary = summary_df,
    effective_sizes = eff_size_df
  )
}
