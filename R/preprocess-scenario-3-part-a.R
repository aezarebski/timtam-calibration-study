suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
library(tidyr)
suppressPackageStartupMessages(library(xml2))
set.seed(1)

## =============================================================================

make_record <- function(x) {
  x |>
    as.numeric() |>
    as.list() |>
    set_names(c("t", "Psi", "X", "Mu", "Omega")) |>
    as.data.frame()
}

process_log_line <- function(l) {
  vs <- l |>
    str_split("\t") |>
    unlist()

  result <- vs[2] |>
    str_split(";") |>
    unlist() |>
    str_replace_all("(t|Psi|X|Mu|Omega)=", "") |>
    str_split(":") |>
    map(make_record) |>
    bind_rows() |>
    mutate(sample = vs[1])

  return(result)
}

## =============================================================================

annotated_tree_and_seqs <- function(recon_tree, sim_data) {
  tmp_tip_labs <- recon_tree$tip.label
  fwd_times_rel <- recon_tree |>
    node.depth.edgelength() |>
    head(length(tmp_tip_labs))

  ix <- nrow(sim_data)
  while (ix > 1) {
    if (sim_data[ix, "Psi"] > sim_data[ix - 1, "Psi"]) {
      last_seq_time <- sim_data[ix, "t"]
      break
    }
    ix <- ix - 1
  }
  rm(ix)

  fwd_times <- last_seq_time + fwd_times_rel - max(fwd_times_rel)
  recon_tree$tip.label <- sprintf("%s_%f", tmp_tip_labs, fwd_times)

  seqs <- simSeq(
    recon_tree,
    l = 1,
    rate = .Machine$double.eps
  )
  return(list(recon_tree = recon_tree, seqs = seqs))
}

## =============================================================================

omega_and_nu_strings <- function(sim_data, sim_duration) {

  as_ssv_str <- function(xs) {
    paste0(as.character(xs), collapse = " ")
  }

  omega_mask <- c(0, diff(sim_data$Omega)) == 1
  omega_times <- sim_duration - sim_data[omega_mask, ]$t
  omega_str <- as_ssv_str(omega_times)
  nu_tbl <- table(floor(sim_duration - sim_data[omega_mask, ]$t))
  nu_ints <- as.integer(nu_tbl)
  nu_times <- as.integer(names(nu_tbl)) + 0.5

  tmp <-
    data.frame(
      nu_times = seq(from = 0, to = sim_duration - 1, by = 1.0) + 0.5,
      nu_default = 0
    ) |>
    left_join(
      data.frame(nu_times = nu_times, nu_ints = nu_ints),
      by = "nu_times"
    ) |>
    mutate(nu_int_full = ifelse(is.na(nu_ints), nu_default, nu_ints))

  list(omega = omega_str,
       nu_times = as_ssv_str(tmp$nu_times),
       nu_counts = as_ssv_str(tmp$nu_int_full))
}

## =============================================================================

remaster_log_file <- "out/s3/remaster-scenario-3.log"
if (!file.exists(remaster_log_file)) {
  stop("The expected remaster log file does not exist")
}

## NOTE that there have been some changes to the output from remaster,
## so the following will try to parse it in one of two ways, hopefully
## one of them will work.
sim_dfs <- tryCatch(
  expr = {
    tail(readLines(remaster_log_file), -1) |>
      map(process_log_line) |>
      bind_rows()},
  error = function(e) {
    warning("It looks like you are using a newer version of remaster. Parsing for remaster.v2.5.1")
    read.table(
      file = remaster_log_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    ) |>
      select(-index) |>
      rename(sample = Sample) |>
      pivot_wider(
        names_from = population,
        values_from = value,
        id_cols = c(sample, t)) |>
      as.data.frame()
  })

num_samples <- sim_dfs |>
  pluck("sample") |>
  unique() |>
  length()

with(
  list(build_node = read_xml("build.xml")),
  {
    num_sims <- build_node |>
      xml_find_first(xpath = "//property[@name='numSims']") |>
      xml_attr("value") |>
      as.integer()
    if (num_sims != num_samples) {
      stop("The number of simulations found in the log file does not match the number requested in build.xml") # nolint
    } else {
      message("The number of simulations found in the log file matches the number requested in build.xml") # nolint
    }
  }
)

write.table(
  x = sim_dfs,
  file = "out/s3/trajectories-scenario-3.csv",
  sep = ",",
  row.names = FALSE
)

recon_trees <- "out/s3/remaster-scenario-3.tree" |>
  readLines() |>
  str_subset("^tree STATE_[0-9]+") |>
  str_replace("tree STATE_[0-9]+ = ", "") |>
  read.tree(text = _)

remaster_node <- "xml/remaster-scenario-3.xml" |> read_xml()
sim_duration <- remaster_node |>
  xml_find_first("//trajectory") |>
  xml_attr("maxTime") |>
  as.numeric()

write.table(
  x = filter(sim_dfs, t == sim_duration),
  file = "out/s3/final-simulation-state.csv",
  sep = ",",
  row.names = FALSE
)

for (ix in seq.int(num_samples)) {
  sim_data <- filter(sim_dfs, sample == as.character(ix - 1))
  tmp <- annotated_tree_and_seqs(recon_trees[[ix]], sim_data)
  om_and_nu <- omega_and_nu_strings(sim_data, sim_duration)

  writeLines(
    text = om_and_nu$omega,
    con = str_interp("out/s3/occurrence-times-scenario-3-sample-$[03d]{ix}.ssv")
  )
  writeLines(
    text = c(om_and_nu$nu_times, om_and_nu$nu_counts),
    con = str_interp(
      "out/s3/nu-times-and-counts-scenario-3-sample-$[03d]{ix}.ssv")
  )

  write.phyDat(
    tmp$seqs,
    str_interp("out/s3/sequences-scenario-3-sample-$[03d]{ix}.fasta"),
    format = "fasta"
  )
  write.tree(
    tmp$recon_tree,
    str_interp("out/s3/reconstruction-scenario-3-sample-$[03d]{ix}.tree")
  )
}
