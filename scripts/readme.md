# Scripts

Run scripts from anywhere with `Rscript path/to/script.R`; each script locates the
project root before sourcing `myfunctions.R`.

- `figures/`: figure-generation scripts. Outputs go to `Figures/`.
- `comparisons/`: benchmark/comparison scripts. Outputs go to `results/`.

Manuscript figure mapping:

- `scripts/figures/Figure01.R` writes `Figures/Figure01.png`.
- `scripts/figures/Figure02.R` writes `Figures/Figure02.png`.
- `scripts/figures/Figure03.R` writes `Figures/Figure03.png`.
- `scripts/figures/Figure04.R` writes `Figures/Figure04.png`.
