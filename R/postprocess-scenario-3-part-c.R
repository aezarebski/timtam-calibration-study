suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(reshape2))
library(xtable)
set.seed(1)

## ================================================================
## The `helper-functions.R` file provides the following functions:
##
## - `read_mcmc`: Read an MCMC log file into a `coda` object.
## - `make_summary`: Summarise an MCMC object with effective size
## - `in_ci`: Check if a value is within the 95% credible interval
## - `ci_width`: Calculate the width of the 95% credible interval
##
source("R/helper-functions.R")
## ================================================================

#' Read a BEAST2 log file into a data frame.
#'
#' @param filename is the path to the log file.
#' @param burn is the number to remove from the start.
#' @param take_last is the number to take from the end.
#'
#' @return data frame containing the samples.
#'
read_beast2_log <- function(filename, burn = 0, take_last = NA) {
  y <- read.csv(filename, sep = "\t", comment.char = "#")
  if (is.na(take_last) && burn >= 0) {
    return(tail(y, nrow(y) - burn))
  } else if (!is.na(take_last) && burn == 0) {
    return(tail(y, take_last))
  } else {
    stop("Unsupported arguments given to read_beast2_log.")
  }
}

remaster_node <- "xml/remaster-scenario-3.xml" |> read_xml()
build_node <- read_xml("build.xml")
num_replicates <- build_node |>
  xml_find_first(xpath = "//property[@name='numSims']") |>
  xml_attr("value") |>
  as.integer()

extract_rate <- function(rate_name) {
  xpath <- stringr::str_interp("//reaction[@id=\"${rate_name}\"]")
  remaster_node |>
    xml_find_first(xpath) |>
    xml_attr("rate") |>
    str_split(" ") |>
    as_vector() |>
    as.numeric()
}

params <- list(
  br1 = extract_rate("lambdaReaction")[1],
  br2 = extract_rate("lambdaReaction")[2],
  dr = extract_rate("muReaction"),
  sr = extract_rate("psiReaction"),
  or = extract_rate("omegaReaction"),
  hs = 0,
  ts_interval = 1.0
)

log_files <-
  list.files("out/s3", pattern = "*.log", full.names = TRUE) |>
  str_subset("3-2")

if (length(log_files) != num_replicates) {
  stop("Number of log files does not match number of replicates.")
}

final_prev_df <- read.csv("out/s3/final-simulation-state.csv") |>
  mutate(replicate = sample + 1) |>
  select(X, replicate)

bias <- function(x, posts) {
  median((posts - x) / x)
}
error <- function(x, posts) {
  median(abs((posts - x) / x))
}
log_file_summary <- function(log_file, params) {
  replicate_num <- log_file |>
    str_extract("[0-9]{3}") |>
    as.integer()
  final_prev <- final_prev_df |>
    filter(replicate == replicate_num) |>
    pluck("X")
  with(read_beast2_log(log_file), {
    net_removal_rate <- params$dr + TTSamplingRate + TTOmegaRate
    TTREff1 <- TTBirthRate.1 / net_removal_rate
    TTREff2 <- TTBirthRate.2 / net_removal_rate
    true_reff_1 <- params$br1 / (params$dr + params$sr + params$or)
    true_reff_2 <- params$br2 / (params$dr + params$sr + params$or)
    data.frame(
      m_birth_rate_1 = median(TTBirthRate.1),
      m_birth_rate_2 = median(TTBirthRate.2),
      m_sampling_rate = median(TTSamplingRate),
      m_occurrence_rate = median(TTOmegaRate),
      m_reff_1 = median(TTREff1),
      m_reff_2 = median(TTREff2),
      b_birth_rate_1 = bias(params$br1, TTBirthRate.1),
      b_birth_rate_2 = bias(params$br2, TTBirthRate.2),
      b_sampling_rate = bias(params$sr, TTSamplingRate),
      b_occurrence_rate = bias(params$or, TTOmegaRate),
      b_reff_1 = bias(true_reff_1, TTREff1),
      b_reff_2 = bias(true_reff_2, TTREff2),
      b_prev = bias(final_prev, HistorySizes),
      e_birth_rate_1 = error(params$br1, TTBirthRate.1),
      e_birth_rate_2 = error(params$br2, TTBirthRate.2),
      e_sampling_rate = error(params$sr, TTSamplingRate),
      e_occurrence_rate = error(params$or, TTOmegaRate),
      e_reff_1 = error(true_reff_1, TTREff1),
      e_reff_2 = error(true_reff_2, TTREff2),
      e_prev = error(final_prev, HistorySizes),
      w_birth_rate_1 = ci_width(params$br1, TTBirthRate.1),
      w_birth_rate_2 = ci_width(params$br2, TTBirthRate.2),
      w_sampling_rate = ci_width(params$sr, TTSamplingRate),
      w_occurrence_rate = ci_width(params$or, TTOmegaRate),
      w_reff_1 = ci_width(true_reff_1, TTREff1),
      w_reff_2 = ci_width(true_reff_2, TTREff2),
      in_ci_birth_rate_1 = in_ci(params$br1, TTBirthRate.1),
      in_ci_birth_rate_2 = in_ci(params$br2, TTBirthRate.2),
      in_ci_sampling_rate = in_ci(params$sr, TTSamplingRate),
      in_ci_occurrence_rate = in_ci(params$or, TTOmegaRate),
      in_ci_reff_1 = in_ci(true_reff_1, TTREff1),
      in_ci_reff_2 = in_ci(true_reff_2, TTREff2),
      in_ci_prev = in_ci(final_prev, HistorySizes),
      log_file = log_file
    )
  })
}

replicate_summaries <- map(log_files, \(lf) log_file_summary(lf, params)) |>
  bind_rows()

result <- with(replicate_summaries, {
  percent_true <- function(bs) {
    100 * sum(bs) / length(bs)
  }
  data.frame(
    parameter = c(
      "birth rate 1", "birth rate 2",
      "sampling rate", "occurrence rate",
      "R effective 1", "R effective 2",
      "prevalence"
    ),
    true = with(params, c(br1, br2, sr, or, br1/(dr+sr+or), br2/(dr+sr+or), hs)),
    median = c(
      median(m_birth_rate_1), median(m_birth_rate_2),
      median(m_sampling_rate), median(m_occurrence_rate),
      median(m_reff_1), median(m_reff_2),
      NA
    ),
    error = c(
      median(e_birth_rate_1), median(e_birth_rate_2),
      median(e_sampling_rate), median(e_occurrence_rate),
      median(e_reff_1), median(e_reff_2),
      median(e_prev)
    ),
    bias = c(
      median(b_birth_rate_1), median(b_birth_rate_2),
      median(b_sampling_rate), median(b_occurrence_rate),
      median(b_reff_1), median(b_reff_2),
      median(b_prev)
    ),
    ci_width = c(
      median(w_birth_rate_1), median(w_birth_rate_2),
      median(w_sampling_rate), median(w_occurrence_rate),
      median(w_reff_1), median(w_reff_2),
      NA
    ),
    ci_percent = c(
      percent_true(in_ci_birth_rate_1), percent_true(in_ci_birth_rate_2),
      percent_true(in_ci_sampling_rate), percent_true(in_ci_occurrence_rate),
      percent_true(in_ci_reff_1), percent_true(in_ci_reff_2),
      percent_true(in_ci_prev)
    )
  )
})
row.names(result) <- NULL

result <- xtable(result)
align(result) <- xalign(result)
digits(result) <- xdigits(result, zap = 3)
display(result) <- xdisplay(result)
print(result,
      include.rownames = FALSE,
      file = "out/s3/summary-scenario-3-2.tex")

#' Hypothesis test on the binomial probability.
#'
#' @param x_obs integer number of successes observed.
#' @param size integer number of tests carried out.
#' @param prob probability of success in each test
#'
calibration_test <- function(x_obs, size, prob = 0.95) {
  pmf_vals <- dbinom(x = 0:size, prob = prob, size = size)
  obs_prob <- dbinom(x = x_obs, prob = prob, size = size)
  mask <- pmf_vals <= obs_prob
  p_val <- sum(pmf_vals[mask])
  list(
    x = x_obs,
    n = size,
    p_val = p_val,
    reject_null = p_val < 0.05
  )
}

list(
  as.data.frame(calibration_test(90, 100)),
  as.data.frame(calibration_test(91, 100)),
  as.data.frame(calibration_test(99, 100)),
  as.data.frame(calibration_test(100, 100))
) |>
  bind_rows() |>
  print()
