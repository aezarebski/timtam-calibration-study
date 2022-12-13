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
  processx::run(
    command = exe, args = args,
    timeout = timeout,
    error_on_status = FALSE,
    cleanup_tree = TRUE
  )
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

simulate_parameters <- function(uid, config) {
  set.seed(uid)
  r0 <- config$param$r0()
  sigma <- config$param$sigma()
  props <- rdirichlet(config$params$removal_weights) # nolint
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

simulate_epidemic <- function(uid, params, config) {
  data_string <-
    stringr::str_interp(
               "simId=$[d]{uid},duration=$[f]{duration},lambda=$[.5f]{lambda},mu=$[.5f]{mu},psi=$[.5f]{psi},omega=$[.5f]{omega},loggerTreeFile=${treefile},loggerLogFile=${logfile}", # nolint
               as.environment(
                 c(params,
                   list(uid = uid,
                        treefile = config$files$simulation$treefile(uid),
                        logfile = config$files$simulation$logfile(uid)))
               )
             )
  ps_return_val <- run_beast( # nolint
    config$files$simulation$remaster_xml,
    config$sim_timeout_seconds,
    c("-D", data_string)
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
