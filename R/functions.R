library(ggplot2)
library(cowplot)

#' Sample from Dirichlet distribution.
#'
#' https://en.wikipedia.org/wiki/Dirichlet_distribution#From_gamma_distribution
#'
rdirichlet <- function(shape) {
  tmp <- rgamma(n = length(shape), shape = shape, rate = 1)
  tmp / sum(tmp)
}

#' Run BEAST2.7.x XML.
run_beast <- function(beast_xml, timeout, ...) {
  exe <- "./lib/beast/bin/beast"
  if (!file.exists(exe)) {
    stop(sprintf("Could not find executable: %s", exe))
  }
  if (!file.exists(beast_xml)) {
    stop(sprintf("Could not find BEAST XML: %s", beast_xml))
  }
  args <- c(..., beast_xml)
  process_val <-
    processx::run(
                command = exe, args = args,
                timeout = timeout,
                error_on_status = FALSE,
                cleanup_tree = TRUE
              )
  return(process_val)
}

#' Read a remaster simulation log file into a sensible data frame.
read_remaster_simulation <- function(logfile) {
  char_to_df <- \(x) {
    eval(
      parse(
        text = paste0(
          c("data.frame(", paste0(tolower(x), collapse = ","), ")"),
          collapse = ""
        )
      )
    )
  }

  logfile |>
    readLines() |>
    purrr::pluck(2) |>
    stringr::str_split("\t") |>
    unlist() |>
    tail(1) |>
    stringr::str_split(";") |>
    unlist() |>
    stringr::str_split(":") |>
    purrr::map(char_to_df) |>
    dplyr::bind_rows()
}

make_uids <- function(num_sims) {
  seq.int(from = 1, to = num_sims)
}

simulate_parameters_a <- function(uid, config) {
  r0 <- config$param$r0()
  sigma <- config$param$sigma()
  if (!any(is.na(config$param$removal_weights))) {
    props <- rdirichlet(config$param$removal_weights)
  } else {
    ## TODO Clean up this brittle code to work with any of the entries being NA.
    tmp_ws <- config$param$removal_weights
    if (all(is.na(tmp_ws) == c(FALSE, FALSE, TRUE))) {
      props <- c(rdirichlet(tmp_ws[c(1, 2)]), 0)
    } else {
      stop(sprintf("Failed to simulate parameters for UID: %d", uid))
    }
  }
  lambda <- r0 * sigma
  mu <- sigma * props[1]
  psi <- sigma * props[2]
  omega <- sigma * props[3]
  list(
    sigma = sigma,
    r0 = r0,
    lambda = lambda,
    mu = mu,
    psi = psi,
    omega = omega,
    seq_length = config$param$seq_len(),
    seq_rate = config$param$seq_rate(),
    duration = config$param$duration()
  )
}

simulate_parameters_c <- function(uid, config) {
  r01 <- config$param$r01()
  r02 <- config$param$r02()
  sigma <- config$param$sigma()
  if (!any(is.na(config$param$removal_weights))) {
    props <- rdirichlet(config$param$removal_weights)
  } else {
    ## TODO Clean up this brittle code to work with any of the entries being NA.
    tmp_ws <- config$param$removal_weights
    if (all(is.na(tmp_ws) == c(FALSE, FALSE, TRUE))) {
      props <- c(rdirichlet(tmp_ws[c(1, 2)]), 0)
    } else {
      stop(sprintf("Failed to simulate parameters for UID: %d", uid))
    }
  }
  lambda1 <- r01 * sigma
  lambda2 <- r02 * sigma
  mu <- sigma * props[1]
  psi <- sigma * props[2]
  omega <- sigma * props[3]
  list(
    sigma = sigma,
    r01 = r01,
    r02 = r02,
    lambda1 = lambda1,
    lambda2 = lambda2,
    lambdaChangeTime = config$param$changeTime(),
    mu = mu,
    psi = psi,
    omega = omega,
    seq_length = config$param$seq_len(),
    seq_rate = config$param$seq_rate(),
    duration = config$param$duration()
  )
}

#' A random set of parameters as specified by the configuration provided. Since
#' we are setting the seed within this function it should be deterministic.
simulate_parameters <- function(uid, config) {
  set.seed(uid)
  remaster_xml <- basename(config$files$simulation$remaster_xml)
  if (remaster_xml == "remaster-a.xml") {
    simulate_parameters_a(uid, config)
  } else if (remaster_xml == "remaster-c.xml") {
    simulate_parameters_c(uid, config)
  }
}

remaster_data_string_a <- function(uid, params, config) {
  stringr::str_interp(
             "simId=$[d]{uid},duration=$[f]{duration},lambda=$[.5f]{lambda},mu=$[.5f]{mu},psi=$[.5f]{psi},omega=$[.5f]{omega},loggerTreeFile=${treefile},loggerLogFile=${logfile}", # nolint
             as.environment(
               c(params,
                 list(uid = uid,
                      treefile = config$files$simulation$treefile(uid),
                      logfile = config$files$simulation$logfile(uid)))
             )
           )
}

remaster_data_string_c <- function(uid, params, config) {
  stringr::str_interp(
             "simId=$[d]{uid},duration=$[f]{duration},lambdaChangeTime=$[f]{lambdaChangeTime},lambda2=$[.5f]{lambda2},lambda1=$[.5f]{lambda1},mu=$[.5f]{mu},psi=$[.5f]{psi},omega=$[.5f]{omega},loggerTreeFile=${treefile},loggerLogFile=${logfile}", # nolint
             as.environment(
               c(params,
                 list(uid = uid,
                      treefile = config$files$simulation$treefile(uid),
                      logfile = config$files$simulation$logfile(uid)))
             )
           )
}

remaster_data_string <- function(uid, params, config) {
  remaster_xml <- basename(config$files$simulation$remaster_xml)
  if (remaster_xml == "remaster-a.xml") {
    remaster_data_string_a(uid, params, config)
  } else if (remaster_xml == "remaster-c.xml") {
    remaster_data_string_c(uid, params, config)
  }
}

simulate_epidemic <- function(uid, params, config) {
  data_string <- remaster_data_string(uid, params, config)
  ps_return_val <- run_beast( # nolint
    config$files$simulation$remaster_xml,
    config$sim_timeout_seconds,
    c("-overwrite", "-D", data_string)
  )
  if (ps_return_val$status == 1 && !ps_return_val$timeout) {
    stop(ps_return_val$stderr)
  }
  if (!ps_return_val$timeout) {
    sim_data <- read_remaster_simulation(config$files$simulation$logfile(uid))
    recon_tree <- ape::read.nexus(config$files$simulation$treefile(uid))
    result <- list(
      simulation = sim_data,
      recon_tree = recon_tree
    )
  } else {
    result <- sprintf("Simulation number %d timed out.", uid)
  }
  return(result)
}

#' Create a simple summary of the epidemic.
summarise_epidemic <- function(uid, epi) {
  if (is.null(epi)) {
    return(epi)
  } else {
    sim_df <- epi$simulation
    list(
      final_prevalence = tail(sim_df$x, 1),
      total_sequences = tail(sim_df$psi, 1),
      total_occurrences = tail(sim_df$omega, 1)
    )
  }
}

#' The multiple sequence alignment on the reconstructed tree.
simulate_genomes <- function(epi, params) {
  if (is.null(epi)) {
    return(epi)
  }

  sim_data <- epi$simulation
  recon_tree <- epi$recon_tree
  tmp_tip_labs <- recon_tree$tip.label
  fwd_times_rel <- recon_tree |>
    ape::node.depth.edgelength() |>
    head(length(tmp_tip_labs))

  ix <- nrow(sim_data)
  while (ix > 1) {
    if (sim_data[ix, "psi"] > sim_data[ix - 1, "psi"]) {
      last_seq_time <- sim_data[ix, "t"]
      break
    }
    ix <- ix - 1
  }
  rm(ix)

  fwd_times <- last_seq_time + fwd_times_rel - max(fwd_times_rel)
  recon_tree$tip.label <- sprintf("%s_%f", tmp_tip_labs, fwd_times)
  phangorn::simSeq(recon_tree, l = params$seq_length, rate = params$seq_rate)
}

#' Compute a daily time series of confirmed cases.
#'
#' Aggregate the point process of unsequenced unscheduled observations into a
#' time series which we can model as scheduled unsequenced observations. Note
#' that the times in the result are backwards and relative to the timing of the
#' final sequenced sample.
extract_time_series <- function(epi, msa) {
  if (is.null(epi)) {
    return(epi)
  }
  if (tail(epi$simulation$omega, 1) < 1) {
    stop("Invalid simulation data! Check the number of omega events")
  }

  occ_df <- data.frame(t = epi$simulation$t, cum_occ = epi$simulation$omega)
  occ_df$inst_occ <- c(0, diff(occ_df$cum_occ))
  ## The following line computes the daily number of cases, this relies on a
  ## resolution of one unit of time!
  agg_occs <- table(ceiling(occ_df[occ_df$inst_occ == 1, ]$t))
  agg_df <- data.frame(
    time = as.integer(names(agg_occs)),
    num_occ = as.integer(agg_occs)
  )

  padding_df <- data.frame(time = seq.int(max(agg_df$time)), empty_obs = 0)
  agg_occ_df <- dplyr::full_join(x = padding_df, y = agg_df, by = "time") |>
    dplyr::mutate(size = ifelse(is.na(num_occ), empty_obs, num_occ)) |> # nolint
    dplyr::select(time, size) # nolint

  last_sample_time <- msa |>
    labels() |>
    stringr::str_match(pattern = "[0-9.]+$") |>
    as.numeric() |>
    max()

  data.frame(
    times = last_sample_time - agg_occ_df$time,
    sizes = agg_occ_df$size
  )
}

#' Turn a numeric vector into a space separated string of values.
numeric_as_ssv <- function(xs) {
  if (is.numeric(xs)) {
    stringr::str_flatten(string = as.character(xs), collapse = " ")
  } else {
    stop("non-numeric value given to numeric_as_ssv.")
  }
}

#' Construct an XML object which could be used by BEAST2 to run an MCMC analysis
#' with timtam.
get_beast_mcmc_xml <- function(uid, msa, ts_df, params, config) {
  beast_root <- xml2::read_xml(config$files$mcmc$template)
  data_node <- xml2::xml_find_first(beast_root, "//data")

  for (ix in seq.int(msa)) {
    leaf_name <- names(msa[ix])
    leaf_seq <- stringr::str_flatten(as.character(msa[ix]))
    seq_node <- xml2::read_xml("<sequence spec='beast.base.evolution.alignment.Sequence' totalcount='4'></sequence>") # nolint
    xml2::xml_set_attr(
      seq_node, attr = "id", value = stringr::str_flatten(c("seq_", leaf_name))
    )
    xml2::xml_set_attr(seq_node, attr = "taxon", value = leaf_name)
    xml2::xml_set_attr(seq_node, attr = "value", value = leaf_seq)
    xml2::xml_add_child(data_node, seq_node)
  }

  if (!is.null(ts_df)) {
    d_time_node <- xml2::xml_find_first(beast_root, "//disasterTimes")
    xml2::xml_set_text(x = d_time_node, value = numeric_as_ssv(ts_df$times))
    d_size_node <- xml2::xml_find_first(beast_root, "//disasterSizes")
    xml2::xml_set_text(x = d_size_node, value = numeric_as_ssv(ts_df$sizes))
  }

  trait_node <- xml2::xml_find_first(beast_root, "//state/tree/trait")
  trait_fn <-
    \(name) sprintf("%s=%s", name, stringr::str_extract(name, "[.0-9]+$"))
  trait_val <-
    stringr::str_flatten(purrr::map_chr(names(msa), trait_fn), collapse = ",")
  xml2::xml_set_attr(trait_node, attr = "value", value = trait_val)

  sigma_node <- xml2::xml_find_first(beast_root, "//distribution/sigma")
  xml2::xml_set_text(sigma_node, value = as.character(params$sigma))

  origin_time_node <-
    xml2::xml_find_first(beast_root, "//distribution/originTime")
  sim_duration <- params$duration
  xml2::xml_set_text(origin_time_node, value = as.character(sim_duration))
  init_tree_node <-
    xml2::xml_find_first(beast_root, "//init[@spec='RandomTree']")
  xml2::xml_set_attr(
    init_tree_node,
    attr = "rootHeight",
    value = as.character(sim_duration - 0.1)
    )
  xml2::xml_find_first(beast_root, "//logger") |>
    xml2::xml_set_attr(attr = "fileName",
                       value = config$files$mcmc$log_file(uid))
  xml2::xml_find_first(beast_root, "//run") |>
    xml2::xml_set_attr(
      attr = "chainLength",
      value = config$mcmc$length
    )
  return(beast_root)
}

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

run_beast_mcmc <- function(beast_xml, config) {
  if (is.element(el = "xml_document", class(beast_xml))) {
    filepath_beast_xml <- tempfile(pattern = "beast-mcmc-", fileext = ".xml")
    xml2::write_xml(x = beast_xml, file = filepath_beast_xml)
  } else if (file.exists(beast_xml)) {
    filepath_beast_xml <- beast_xml
  } else {
    stop("run_beast_mcmc could not work out how to handle beast_xml argument.")
  }
  run_beast(filepath_beast_xml, config$mcmc_timeout_seconds, c("-overwrite"))

  log_file <-
    xml2::xml_attr(xml2::xml_find_first(beast_xml, "//logger"), "fileName")
  read_beast2_log(log_file, burn = config$mcmc$num_burn) |>
    dplyr::mutate(log_file = log_file) |>
    dplyr::select(-Sample) # nolint
}

est_and_cri <- function(x, p) {
  tail_weight <- 0.5 * (1 - p)
  quantile(x, probs = c(tail_weight, 0.5, 1 - tail_weight))
}

summarise_post_samples <- function(uid, mcmc_samples, epi_summary, params,
                                   config) {
  coda_mcmc <-
    coda::as.mcmc(dplyr::select(mcmc_samples, -log_file)) # nolint
  r0_ci <- est_and_cri(mcmc_samples$TTR0, config$cri_coverage)
  final_prev <- epi_summary$final_prevalence
  prev_ci <- est_and_cri(mcmc_samples$TTHistorySizes, config$cri_coverage)

  data.frame(
    uid = uid,
    min_ess = min(coda::effectiveSize(coda_mcmc)),
    num_samples = nrow(mcmc_samples),
    ci_prob = config$cri_coverage,
    r0_true = params$r0,
    r0_lower = r0_ci[1],
    r0_point = r0_ci[2],
    r0_upper = r0_ci[3],
    contains_r0 = r0_ci[1] < params$r0 & params$r0 < r0_ci[3],
    prevalence_true = final_prev,
    prevalence_lower = prev_ci[1],
    prevalence_point = prev_ci[2],
    prevalence_upper = prev_ci[3],
    contains_prevalence = (prev_ci[1] < final_prev & final_prev < prev_ci[3])
  ) |>
    dplyr::bind_cols(as.data.frame(epi_summary))
}

#' Plot the prevalence estimates against their true values.
plot_prevalence_estimates <- function(mcmc_summary) {
  ggplot() +
    geom_abline(intercept = 0, slope = 1) +
    geom_linerange(
      data = mcmc_summary,
      mapping = aes(
        x = prevalence_true, # nolint
        ymin = prevalence_lower, # nolint
        ymax = prevalence_upper # nolint
      ),
      linewidth = 2,
      alpha = 0.3
    ) +
    geom_point(
      data = mcmc_summary,
      mapping = aes(
        x = prevalence_true,
        y = prevalence_point, # nolint
        shape = contains_prevalence # nolint
      )
    ) +
    scale_y_log10() +
    scale_x_log10() +
    scale_shape_manual(values = c(4, 16), breaks = c(FALSE, TRUE)) +
    labs(x = "Prevalence", y = "Estimated prevalence") +
    theme_bw()
}

#' Plot the R-naught estimates against their true values.
plot_r0_estimates <- function(mcmc_summary) {
  ggplot() +
    geom_abline(intercept = 0, slope = 1) +
    geom_linerange(
      data = mcmc_summary,
      mapping = aes(
        x = r0_true, # nolint
        ymin = r0_lower, # nolint
        ymax = r0_upper  # nolint
      ),
      linewidth = 2,
      alpha = 0.3
    ) +
    geom_point(
      data = mcmc_summary,
      mapping = aes(
        x = r0_true,
        y = r0_point, # nolint
        shape = contains_r0 # nolint
      ),
      size = 2
    ) +
    scale_shape_manual(values = c(4, 16), breaks = c(FALSE, TRUE)) +
    labs(x = "R-naught", y = "Estimated R-naught") +
    theme_bw()
}

#' Write a ggplot to a file in an A6 size.
save_plot_a6 <- function(gg, fn) {
 ggsave(filename = fn,
        plot = gg,
        height = 10.5, width = 14.8,
        units = "cm")
 return(fn)
}

combined_figure <- function(prev_gg, r0_gg) {
  legend <-
    get_legend(
               r0_gg +
               labs(shape = "CrI contains true value") +
               theme(legend.position = "top",
                     legend.box.margin = margin(0, 0, 0, 12))
             )

  plot_grid(
    legend,
    prev_gg + theme(legend.position = "none"),
    r0_gg + theme(legend.position = "none"),
    ncol = 1,
    rel_heights = c(0.25, 1, 1),
    labels = c("", "A", "B")
  )
}
