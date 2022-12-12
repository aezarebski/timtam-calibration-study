#' Run BEAST2.7.x XML.
run_beast <- function(beast_xml, timeout, ...) {
  exe <- "./lib/beast/bin/beast"
  if (!file.exists(exe)) {
    stop(sprintf("Could not find executable: %s", exe))
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
               "simId=$[d]{uid},duration=$[f]{duration},lambda=$[.5f]{lambda},mu=$[.5f]{mu},psi=$[.5f]{psi},omega=$[.5f]{omega}", # nolint
               as.environment(c(params, list(uid = uid)))
             )
  ps_return_val <- run_beast( # nolint
    config$files$simulation$remaster_xml,
    config$sim_timeout_seconds,
    c("-D", data_string)
  )
  if (ps_return_val$status == 1 & !ps_return_val$timeout) {
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
