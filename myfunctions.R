.exit_bessel_root <- local({
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
    dirname(normalizePath(ofiles[[length(ofiles)]], mustWork = TRUE))
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
})

.exit_bessel_files <- file.path(
  .exit_bessel_root,
  "R",
  c(
    "validation.R",
    "distributions.R",
    "bounds_g.R",
    "sampler_g.R",
    "bounds_f.R",
    "envelope_f.R",
    "sampler_f.R"
  )
)

.exit_bessel_env <- environment()

for (.exit_bessel_file in .exit_bessel_files) {
  source(.exit_bessel_file, local = .exit_bessel_env)
}

rm(.exit_bessel_file, .exit_bessel_files, .exit_bessel_env, .exit_bessel_root)
