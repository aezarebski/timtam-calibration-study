## This R script provides some helpful functions which are used across
## other scripts. However, substantial care must be taken when
## sourcing this file as it introduces multiple functions so you need
## to make sure it will not overwrite any existing functions in your
## environment. For this reason it is a good idea to source it at the
## start of your script and include the following comments to make it
## clear what is happening:
##
## ## ================================================================
## ## The `helper-functions.R` file provides the following functions:
## ##
## ## - `read_mcmc`: Read an MCMC log file into a `coda` object.
## ## - `make_summary`: Summarise an MCMC object with effective size
## ## - `in_ci`: Check if a value is within the 95% credible interval
## ## - `ci_width`: Calculate the width of the 95% credible interval
## ## - `my_ggsave`: Save a ggplot as PNG and SVG in one call
## ##
## source("R/helper-functions.R")
## ## ================================================================
##
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(ggplot2))


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
  ## summary_df <- summary(mcmc_obj)$quantiles
  ## tmp <- rownames(summary_df)
  ## summary_df <- summary_df |>
  ##   as.data.frame() |>
  ##   mutate(variable = tmp, replicate = ix) |>
  ##   rename(
  ##     fns_1 = "2.5%",
  ##     fns_2 = "25%",
  ##     fns_3 = "50%",
  ##     fns_4 = "75%",
  ##     fns_5 = "97.5%",
  ##   )
  cri_95 <- HPDinterval(mcmc_obj, prob = 0.95)
  cri_50 <- HPDinterval(mcmc_obj, prob = 0.5)
  cri_median <- as.numeric(apply(X = mcmc_obj, FUN = median, MARGIN = 2))
  summary_df <- data.frame(
    fns_1 = cri_95[, "lower"],
    fns_2 = cri_50[, "lower"],
    fns_3 = cri_median,
    fns_4 = cri_50[, "upper"],
    fns_5 = cri_95[, "upper"],
    variable = colnames(mcmc_obj),
    replicate = ix
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


#' Check if credible interval contains true value
#'
#' Check if true value is within the 95% credible interval of MCMC
#' samples
#'
#' @param x True value
#' @param posts MCMC samples
#' @return Logical
#'
in_ci <- function(x, posts) {
  ## ci <- quantile(posts, probs = c(0.025, 0.975))
  ## ci[1] <= x && x <= ci[2]
  tmp <- bayestestR::hdi(posts, ci = 0.95)
  tmp$CI_low <= x && x <= tmp$CI_high
}


#' Calculate width of 95% credible interval
#'
#' Calculate the width of the 95% credible interval of MCMC samples
#' normalised by the true value
#'
#' @param x True value
#' @param posts MCMC samples
#' @return Width of credible interval
#'
ci_width <- function(x, posts) {
  ## ci <- quantile(posts, probs = c(0.025, 0.975))
  ## (ci[2] - ci[1]) / x
  tmp <- bayestestR::hdi(posts, ci = 0.95)
  (tmp$CI_high - tmp$CI_low) / x
}


#' Save a plot as both a PNG and an SVG file.
#'
#' This function saves a plot as both a PNG and an SVG file. The
#' filename should end with `.png` and the plot will be saved as a PNG
#' file with the same name. The SVG file will have the same name but
#' with `.svg` instead of `.png`.
#'
#' @param filename The name of the file to save the plot to.
#' @param plot The plot to save.
#' @param ... Additional arguments to pass to `ggsave`.
#'
my_ggsave <- function(filename, plot, ...) {
  ggplot2::ggsave(filename = filename, plot = plot, ...)
  if (str_detect(filename, ".png$")) {
    svg_filename <- str_replace(filename, ".png$", ".svg")
    ggplot2::ggsave(filename = svg_filename, plot = plot, ...)
  }
}
