# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**survtrans** — an R package for transfer learning in survival analysis. Transfers survival information from source domain(s) to a target domain using Cox proportional hazards models with global and local shrinkage (ADMM optimization). Supports lasso, SCAD, and MCP penalties.

## Common Commands

```bash
# Build and check
R CMD build .
R CMD check survtrans_*.tar.gz

# Install locally
R CMD INSTALL .

# Run tests (testthat v3)
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-coxtrans.R")'

# Rebuild C++ bindings after modifying src/
Rscript -e 'Rcpp::compileAttributes()'

# Regenerate documentation after modifying roxygen comments
Rscript -e 'devtools::document()'

# Full check (equivalent to CI)
Rscript -e 'devtools::check()'

# Build README from Rmd
Rscript -e 'devtools::build_readme()'

# Lint and format
Rscript -e 'lintr::lint_package()'
Rscript -e 'styler::style_pkg()'

# Pre-commit hooks
pre-commit install
pre-commit run --all-files
```

## Architecture

### Core Functions

- **`coxtrans()`** (`R/coxtrans.R`): Main entry point. Fits a transfer learning Cox model with three penalty parameters: `lambda1` (sparsity), `lambda2` (global bias shrinkage), `lambda3` (local bias shrinkage). Uses ADMM internally. Returns a `coxtrans` S3 object with methods: `coef`, `predict`, `summary`, `print`, `basehaz`, `vcov`, `logLik`, `diagnose`.

- **`ncvcox()`** (`R/ncvcox.R`): Non-convex penalized Cox model (standalone, no transfer). Returns an `ncvcox` S3 object with similar methods.

- **`survtrans_control()`** (`R/survtrans_control.R`): Shared control parameters for the fitting algorithm (max iterations, tolerance, etc.).

### C++ Layer (src/)

Performance-critical computations via Rcpp/RcppArmadillo with OpenMP support:

- `approx_likelihood.cpp` — approximate partial likelihood, weights, residuals
- `calc_grad_hess.cpp` — gradient and Hessian
- `threshold_prox.cpp` — proximal operator for penalty thresholding
- `close_update.cpp` — closed-form ADMM updates

After modifying any `.cpp` file, run `Rcpp::compileAttributes()` to regenerate `RcppExports.cpp` and `R/RcppExports.R`.

### Coefficients Structure

`coxtrans` returns a coefficients matrix with columns `[group1, group2, ..., groupK, Center]`. Each group's full beta = `coefficients[, group] + coefficients[, "Center"]`. The target group is always placed first in the column order.

### Data Flow

1. User calls `coxtrans(formula, data, group, ...)` where `group` identifies source/target domains
2. `preprocess()` in `R/utils.R` handles formula parsing, model matrix construction (`scale(x)`), and group splitting (sorted by time descending)
3. ADMM iterations call C++ functions for likelihood/gradient/proximal steps
4. Post-processing: coefficients are un-standardized via `sweep(coefficients, 1, x_scale, "/")`; stored `x` is centered but not scaled

### Key Utilities (`R/utils.R`)

- `basehaz()` — Breslow baseline hazard estimator
- `simsurv_tl()` — simulate multi-source survival data for experiments
- `calc_lambda_max()` — compute maximum meaningful lambda for penalty path
- `diagnose()` — diagnostic plots for model assessment

### Scripts (`scripts/`)

Experiment scripts for model evaluation (not part of the R package):

- `01_cv_grid_search.R` — 3D grid search CV for coxtrans (lambda1/lambda2/lambda3), uses partial log-likelihood deviance as metric
- `02_compare_ncvreg_coxtrans.R` — compare ncvreg (target only) vs coxtrans (lambda2=lambda3=0) variable selection
