suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xml2))
set.seed(1)


make_record <- function(x) {
  x |>
    as.numeric() |>
    as.list() |>
    set_names(c("t", "Psi", "X", "Mu")) |>
    as.data.frame()
}

process_log_line <- function(l) {
  vs <- l |>
    str_split("\t") |>
    unlist()

  result <- vs[2] |>
    str_split(";") |>
    unlist() |>
    str_replace_all("(t|Psi|X|Mu)=", "") |>
    str_split(":") |>
    map(make_record) |>
    bind_rows() |>
    mutate(sample = vs[1])

  return(result)
}

sim_log_lines <- tail(readLines("remaster-scenario-1.log"), -1)
sim_dfs <- sim_log_lines |> map(process_log_line) |> bind_rows()

write.table(x = sim_dfs,
            file = "trajectories-scenario-1.csv",
            sep = ",",
            row.names = FALSE)

stop()

## TODO Fix this so it will extract the trees and write the results to a
## sensible output file.

recon_tree <- "remaster-scenario-1.tree" |> read.nexus() |> pluck(1)
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
write.phyDat(seqs, "demo-sequences-scenario-1-sample-1.fasta", format = "fasta")
