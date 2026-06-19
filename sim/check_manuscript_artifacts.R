#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
tex_file <- if (length(args) >= 1) args[[1]] else file.path("tex", "BMA.tex")
out_file <- if (length(args) >= 2) args[[2]] else file.path("sim", "output", "tables", "table_manuscript_artifact_check.csv")

if (!file.exists(tex_file)) {
  stop("TeX file not found: ", tex_file, call. = FALSE)
}

lines <- readLines(tex_file, warn = FALSE)

extract_matches <- function(pattern, kind) {
  pieces <- regmatches(lines, gregexpr(pattern, lines, perl = TRUE))
  vals <- unlist(pieces, use.names = FALSE)
  if (!length(vals)) {
    return(data.frame(kind = character(), manuscript_path = character()))
  }
  paths <- sub(pattern, "\\1", vals, perl = TRUE)
  data.frame(kind = kind, manuscript_path = paths, stringsAsFactors = FALSE)
}

inputs <- extract_matches("\\\\input\\{([^}]+)\\}", "table")
graphics <- extract_matches("\\\\includegraphics(?:\\[[^]]*\\])?\\{([^}]+)\\}", "figure")
artifacts <- unique(rbind(inputs, graphics))

check_one <- function(path) {
  primary <- file.path(dirname(tex_file), path)
  root <- path
  primary_exists <- file.exists(primary)
  root_exists <- file.exists(root)
  found <- if (primary_exists) primary else if (root_exists) root else NA_character_
  data.frame(
    primary_path = primary,
    root_path = root,
    exists_primary = primary_exists,
    exists_root = root_exists,
    found_path = found,
    status = if (is.na(found)) "missing" else "found",
    stringsAsFactors = FALSE
  )
}

checked <- do.call(rbind, lapply(artifacts$manuscript_path, check_one))
report <- cbind(
  data.frame(
    checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    tex_file = tex_file,
    stringsAsFactors = FALSE
  ),
  artifacts,
  checked
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(report, out_file, row.names = FALSE)

missing <- report[report$status == "missing", , drop = FALSE]
message("Checked ", nrow(report), " manuscript table/figure artifacts.")
message("Wrote ", out_file)

if (nrow(missing) > 0) {
  message("Missing artifacts:")
  for (i in seq_len(nrow(missing))) {
    message("  - ", missing$manuscript_path[[i]])
  }
  quit(status = 1)
}

message("All manuscript table/figure artifacts were found.")
