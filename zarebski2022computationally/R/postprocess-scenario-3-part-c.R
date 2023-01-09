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
num_replicates <- remaster_node |>
  xml_find_first("//run") |>
  xml_attr("nSims") |>
  as.integer()

log_files <-
  list.files("out/s3", pattern = "*.log", full.names = TRUE) |>
  str_subset("3-2")

if (length(log_files) != num_replicates) {
  stop("Number of log files does not match number of replicates.")
}

log_file_summary <- function(log_file, params) {
  in_ci <- function(x, posts) {
    ci <- quantile(posts, probs = c(0.025, 0.975))
    ci[1] <= x && x <= ci[2]
  }
  ci_width <- function(posts) {
    ci <- quantile(posts, probs = c(0.025, 0.975))
    ci[2] - ci[1]
  }
  bias <- function(x, posts) {
    median((posts - x) / x)
  }
  error <- function(x, posts) {
    median(abs((posts - x) / x))
  }
  with(read_beast2_log(log_file), {
    data.frame(
      m_birth_rate_1 = median(TTBirthRate.1),
      m_birth_rate_2 = median(TTBirthRate.2),
      m_sampling_rate = median(TTSamplingRate),
      m_occurrence_rate = median(TTOmegaRate),
      b_birth_rate_1 = bias(params$br1, TTBirthRate.1),
      b_birth_rate_2 = bias(params$br2, TTBirthRate.2),
      b_sampling_rate = bias(params$sr, TTSamplingRate),
      b_occurrence_rate = bias(params$or, TTOmegaRate),
      e_birth_rate_1 = error(params$br1, TTBirthRate.1),
      e_birth_rate_2 = error(params$br2, TTBirthRate.2),
      e_sampling_rate = error(params$sr, TTSamplingRate),
      e_occurrence_rate = error(params$or, TTOmegaRate),
      w_birth_rate_1 = ci_width(TTBirthRate.1),
      w_birth_rate_2 = ci_width(TTBirthRate.2),
      w_sampling_rate = ci_width(TTSamplingRate),
      w_occurrence_rate = ci_width(TTOmegaRate),
      in_ci_birth_rate_1 = in_ci(params$br1, TTBirthRate.1),
      in_ci_birth_rate_2 = in_ci(params$br2, TTBirthRate.2),
      in_ci_sampling_rate = in_ci(params$sr, TTSamplingRate),
      in_ci_occurrence_rate = in_ci(params$or, TTOmegaRate),
      log_file = log_file
    )
  })
}

params <- list(br1 = 0.185, br2 = 0.0925, sr = 0.008, or = 0.046)
replicate_summaries <- map(log_files, \(lf) log_file_summary(lf, params)) |> bind_rows()

result <- with(replicate_summaries, {
  percent_true <- function(bs) {
    100 * sum(bs) / length(bs)
  }
  data.frame(
    parameter = c(
      "birth rate 1", "birth rate 2",
      "sampling rate", "occurrence rate"
    ),
    true = as_vector(params),
    median = c(
      median(m_birth_rate_1), median(m_birth_rate_2),
      median(m_sampling_rate), median(m_occurrence_rate)
    ),
    error = c(
      median(e_birth_rate_1), median(e_birth_rate_2),
      median(e_sampling_rate), median(e_occurrence_rate)
    ),
    bias = c(
      median(b_birth_rate_1), median(b_birth_rate_2),
      median(b_sampling_rate), median(b_occurrence_rate)
    ),
    ci_width = c(
      median(w_birth_rate_1), median(w_birth_rate_2),
      median(w_sampling_rate), median(w_occurrence_rate)
    ),
    ci_percent = c(
      percent_true(in_ci_birth_rate_1), percent_true(in_ci_birth_rate_2),
      percent_true(in_ci_sampling_rate), percent_true(in_ci_occurrence_rate)
    )
  )
})
row.names(result) <- NULL

result <- xtable(result)
align(result) <- xalign(result)
digits(result) <- xdigits(result, zap = 3)
display(result) <- xdisplay(result)
print(result, include.rownames = FALSE)

## % latex table generated in R 4.2.1 by xtable 1.8-4 package
## % Mon Jan  9 13:39:48 2023
## \begin{table}[ht]
## \centering
## \begin{tabular}{lrrrrrr}
##   \hline
## parameter & true & median & error & bias & ci\_width & ci\_percent \\
##   \hline
## birth rate 1 & 0.1850 & 0.1851 & 0.117 & 0.0004 & 0.0945 & 94 \\
##   birth rate 2 & 0.0925 & 0.0948 & 0.334 & 0.0246 & 0.1273 & 96 \\
##   sampling rate & 0.0080 & 0.0103 & 0.339 & 0.2891 & 0.0145 & 92 \\
##   occurrence rate & 0.0460 & 0.0525 & 0.248 & 0.1408 & 0.0557 & 97 \\
##    \hline
## \end{tabular}
## \end{table}


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
) |> bind_rows() |> print()
