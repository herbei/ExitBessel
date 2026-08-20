exit_bessel_script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  }

  ofiles <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) {
        NA_character_
      } else {
        frame$ofile
      }
    },
    character(1L)
  )
  ofiles <- ofiles[!is.na(ofiles)]

  if (length(ofiles) > 0L) {
    return(normalizePath(ofiles[[length(ofiles)]], mustWork = TRUE))
  }

  stop("Could not determine the running script path.", call. = FALSE)
}

# -----------------------------------------------------------------------------

exit_bessel_project_root <- function(script_file = exit_bessel_script_file()) {
  start <- dirname(normalizePath(script_file, mustWork = TRUE))
  candidates <- unique(c(
    normalizePath(getwd(), mustWork = TRUE),
    normalizePath(start, mustWork = TRUE),
    normalizePath(file.path(start, ".."), mustWork = TRUE),
    normalizePath(file.path(start, "..", ".."), mustWork = TRUE),
    normalizePath(file.path(start, "..", "..", ".."), mustWork = TRUE)
  ))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "myfunctions.R")) &&
        dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
  }

  stop("Could not locate the Exit_Bessel_Rcode project root.", call. = FALSE)
}

# -----------------------------------------------------------------------------

exit_bessel_path <- function(..., root = exit_bessel_project_root()) {
  file.path(root, ...)
}

# -----------------------------------------------------------------------------

source_exit_bessel <- function(root = exit_bessel_project_root(), local = parent.frame()) {
  source(file.path(root, "myfunctions.R"), local = local)
}
