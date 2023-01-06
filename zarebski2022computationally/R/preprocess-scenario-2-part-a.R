suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
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


sim_log_lines <- tail(readLines("out/s2/remaster-scenario-2.log"), -1)
sim_dfs <- sim_log_lines |>
  map(process_log_line) |>
  bind_rows()

num_samples <- sim_dfs |>
  pluck("sample") |>
  unique() |>
  length()

write.table(
  x = sim_dfs,
  file = "out/s2/trajectories-scenario-2.csv",
  sep = ",",
  row.names = FALSE
)

recon_trees <- "out/s2/remaster-scenario-2.tree" |>
  readLines() |>
  str_subset("^tree STATE_[0-9]+") |>
  str_replace("tree STATE_[0-9]+ = ", "") |>
  read.tree(text = _)

remaster_node <- "xml/remaster-scenario-2.xml" |> read_xml()
sim_duration <- remaster_node |>
  xml_find_first("//trajectory") |>
  xml_attr("maxTime") |>
  as.numeric()

for (ix in seq.int(num_samples)) {
  sim_data <- filter(sim_dfs, sample == as.character(ix - 1))
  tmp <- annotated_tree_and_seqs(recon_trees[[ix]], sim_data)

  omega_mask <- c(0, diff(sim_data$Omega)) == 1
  omega_str <- paste0(
    as.character(sim_duration - sim_data[omega_mask, ]$t),
    collapse = " "
  )
  writeLines(
    text = omega_str,
    con = str_interp("out/s2/occurrence-times-scenario-2-sample-$[03d]{ix}.ssv")
  )

  write.phyDat(
    tmp$seqs,
    str_interp("out/s2/sequences-scenario-2-sample-$[03d]{ix}.fasta"),
    format = "fasta"
  )
  write.tree(
    tmp$recon_tree,
    str_interp("out/s2/reconstruction-scenario-2-sample-$[03d]{ix}.tree")
  )
}
