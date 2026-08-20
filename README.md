# Exit Bessel R Code

This folder contains the R code and generated outputs for the manuscript.

## Contents

- `R/`: core sampler, density, bound, validation, and SH26 comparison functions.
- `myfunctions.R`: loader for the core implementation.
- `scripts/figures/`: scripts for manuscript Figures 01--04.
- `scripts/comparisons/`: runtime comparison scripts.
- `results/`: CSV result files used by the runtime comparison figure.
- `Figures/`: current generated manuscript figure PNGs.
- `R/SH26_functions.R`: legacy SH26 helper needed by the comparison scripts.

Excluded material: exploratory checks, archived figure variants, old legacy
figures, R history files, cached RData files, PDFs, and LaTeX build artifacts.

## Source Check

Parse all R source files from this directory:

```sh
Rscript -e 'files <- c(list.files("R", "[.]R$", full.names=TRUE), "myfunctions.R", list.files("scripts", "[.]R$", recursive=TRUE, full.names=TRUE)); invisible(lapply(files, parse)); cat("Parsed", length(files), "R files successfully.\n")'
```

## Regenerating Figures

Run from this directory:

```sh
Rscript scripts/figures/Figure01.R
Rscript scripts/figures/Figure02.R
EXIT_BESSEL_RUNTIME_USE_CACHE=true Rscript scripts/figures/Figure03.R
Rscript scripts/figures/Figure04.R
```

The runtime comparison scripts can be rerun with:

```sh
Rscript scripts/comparisons/comparison_f.R
Rscript scripts/comparisons/comparison_g.R
```

`comparison_f.R` uses `EXIT_BESSEL_B` to control accepted draws per setting;
the default is `10000`.
