# Project instructions

This repository is being converted from an existing collection of R code into a well-structured R package.

The primary goal is to preserve the behavior and scientific/statistical logic of the existing R code while reorganizing it into a maintainable package.

## General approach

- Inspect and understand the existing R code before restructuring it.
- Preserve existing calculations and behavior unless I explicitly request a change.
- Refactor incrementally rather than rewriting everything at once.
- Separate reusable functions from analysis scripts and examples.
- Avoid unnecessary changes that are unrelated to packaging the code.

## R package conventions

Follow standard R package conventions.

Use:
- `DESCRIPTION` for package metadata and dependencies.
- `R/` for package functions.
- `tests/testthat/` for automated tests.
- `man/` for generated documentation.
- `vignettes/` for longer examples or workflows when appropriate.
- `README.md` for package-level usage examples.

Use `roxygen2` comments to document exported functions.

Do not manually edit `NAMESPACE` when it can be generated from roxygen documentation.

Use `Imports`, `Suggests`, and other dependency fields appropriately in `DESCRIPTION`. Do not add dependencies without a clear reason.

Prefer `testthat` for package tests.

## Validation

After meaningful package changes, run appropriate checks when the required R tooling is available.

Prefer, as applicable:
- `devtools::document()`
- `devtools::test()`
- `devtools::check()`

Treat `R CMD check` issues seriously. Do not suppress warnings or notes merely to make checks pass without understanding their cause.

## Refactoring existing code

When converting existing scripts into package functions:

1. Identify reusable logical units.
2. Extract them into functions with explicit arguments rather than relying on global variables.
3. Avoid hidden dependencies on the user's global environment or working directory.
4. Avoid `setwd()` inside package code.
5. Use package-safe methods for accessing files and resources.
6. Preserve numerical/statistical behavior.
7. Add tests that compare the refactored function's output with the existing implementation whenever practical.

## Changes

Before making a large structural change, explain the proposed package structure and why.

When modifying existing functionality, summarize:
- what changed,
- why it changed,
- whether behavior changed,
- and what tests or checks were performed.