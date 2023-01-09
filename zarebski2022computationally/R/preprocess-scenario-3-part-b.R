suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xml2))
set.seed(1)

make_timtam <- function(recon_tree_file, occurrence_txt_file, fasta_file, input_xml) {
  timtam <- read_xml(input_xml)
  data_node <- xml_find_first(timtam, xpath = "//data")
  trait_node <- xml_find_first(timtam, xpath = "//trait[@id=\"dateTrait.t\"]")
  tree_node <- xml_find_first(timtam, xpath = "//init[@id=\"NewickTree.t\"]")
  occ_node <- xml_find_first(timtam, xpath = "//parameter[@id=\"OccurrenceTimes\"]")

  fasta_lines <- readLines(fasta_file)
  taxon_strs <- fasta_lines[seq.int(from = 1, to = length(fasta_lines), by = 2)] |> str_match("[0-9]+_[.0-9]+")
  seq_node_strs <- sprintf("<sequence id=\"seq_%s\" spec=\"Sequence\" taxon=\"%s\" totalcount=\"4\" value=\"t\"/>", taxon_strs, taxon_strs)

  for (jx in seq_along(seq_node_strs)) {
    tmp_node <- read_xml(x = seq_node_strs[jx])
    xml_add_child(data_node, tmp_node)
  }

  trait_str <- taxon_strs |>
    str_replace("([0-9]+)_([.0-9]+$)", "\\1_\\2=\\2") |>
    str_flatten(collapse = ",")
  xml_set_attr(trait_node, attr = "value", value = trait_str)

  newick_str <- recon_tree_file |>
    read.tree() |>
    write.tree()
  xml_set_attr(tree_node, attr = "newick", value = newick_str)

  occ_times_str <- readLines(occurrence_txt_file)
  xml_set_text(occ_node, occ_times_str)

  return(timtam)
}

remaster_node <- "xml/remaster-scenario-3.xml" |> read_xml()
num_replicates <- remaster_node |>
  xml_find_first("//run") |>
  xml_attr("nSims") |>
  as.integer()

for (ix in seq.int(num_replicates)) {
  recon_tree_file <- str_interp("out/s3/reconstruction-scenario-3-sample-$[03d]{ix}.tree")
  occurrence_txt_file <- str_interp("out/s3/occurrence-times-scenario-3-sample-$[03d]{ix}.ssv")
  fasta_file <- str_interp("out/s3/sequences-scenario-3-sample-$[03d]{ix}.fasta")
  output_xml <- str_interp("out/s3/timtam-scenario-3-sample-$[03d]{ix}.xml")

  write_xml(
    make_timtam(
      recon_tree_file,
      occurrence_txt_file,
      fasta_file,
      "xml/timtam-scenario-3-1-empty.xml"
    ),
    str_interp("out/s3/timtam-scenario-3-1-sample-$[03d]{ix}.xml")
  )
  write_xml(
    make_timtam(
      recon_tree_file,
      occurrence_txt_file,
      fasta_file,
      "xml/timtam-scenario-3-2-empty.xml"
    ),
    str_interp("out/s3/timtam-scenario-3-2-sample-$[03d]{ix}.xml")
  )
}
