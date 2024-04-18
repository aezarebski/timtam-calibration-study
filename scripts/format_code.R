#!/usr/bin/env Rscript

library(styler)


#' Style the R files in the specified directory
#'
#' @param dir_path The path to the directory containing the R files
#'
style_r_files <- function(dir_path) {
  if (!grepl("/$", dir_path)) {
    dir_path <- paste0(dir_path, "/")
  }

  files <- list.files(dir_path, pattern = "\\.R$", full.names = TRUE)
  lapply(files, style_file)

  if (length(files) > 0) {
    cat("Styled files:\n", paste(files, collapse = "\n"), "\n")
  } else {
    cat("No R files found in the specified directory.\n")
  }
}


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("No directory specified. Usage: ./stylr.R <path_to_directory>", call. = FALSE)
  }

  style_r_files(args[1])
}

main()
