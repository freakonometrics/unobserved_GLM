# Audit-Certified Fair GLMs — numerical reproduction

This repository contains the code needed to reproduce the simulated and real
data experiments, tables, and figures for **Audit-Certified Fair GLMs with
Partially Observed Demographics**.

The repository is intentionally compact:

```text
reproduce.qmd
R/
  01_core.R
  02_data.R
  03_experiments.R
data/
results/
figures/
simulation_outputs_sameX/
```

The Quarto notebook is the single entry point and renders to GitHub-flavoured
Markdown (`reproduce.md`), not HTML.

## Requirements

Install the R packages:

```r
install.packages(c("readxl", "ggplot2"))
```

Quarto is needed only to render the notebook.

## Quick reproducibility check

```bash
quarto render reproduce.qmd
```

This uses reduced Monte Carlo sizes but executes both the original
source/target proxy-provenance simulation and the UCI Credit application.

## Full main-paper reproduction

```bash
quarto render reproduce.qmd -P quick:false
```

This restores the manuscript Monte Carlo sizes, including:

- proxy provenance: 120,000 source observations, 30,000 target observations,
  500 repeated audits per scenario/rate cell;
- Credit fixed-score estimation: 300 audit draws;
- Credit fair learning: 40 repetitions;
- audit rates 1%, 2%, 5%, 10%, and 20%;
- the 31-point positive lambda grid plus zero.

## Full appendix robustness suite

```bash
quarto render reproduce.qmd \
  -P quick:false \
  -P run_robustness:true
```

The robustness suite covers proxy stress tests, exact versus convex fairness
penalties, analytic variance validation, external-pilot Neyman allocation,
counted-budget two-phase Neyman allocation, and lambda-selection robustness.

## Data

The source/target proxy-provenance experiment is fully synthetic.

The real-data application uses the public UCI **Default of Credit Card
Clients** dataset. `R/02_data.R` downloads the original archive automatically
and performs all preprocessing from the raw `.xls` file.

## Main generated figures

The source/target simulation generates:

```text
figure_sameX_provenance_rmse.pdf
figure_sameX_external_vs_target_ratio.pdf
figure_sameX_external_proxy_bias.pdf
figure_sameX_proxy_brier.pdf
figure_sameX_calibration_variance.pdf
figure_sameX_calibration_plugin_bias.pdf
```

The Credit application generates:

```text
credit_moment_rmse_sex.pdf
credit_selected_dp_sex.pdf
credit_frontier_sex.pdf
```

The optional robustness suite generates the v5.1 appendix figures in the same
`figures/` directory.

## Reproducibility note

All Monte Carlo seeds are fixed in the code. The execution environment used to
assemble this repository did not contain R or Quarto, so the final repository
was checked structurally but not executed here. Running the quick render is the
recommended first verification step.
