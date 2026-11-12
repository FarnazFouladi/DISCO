# Rebuilding the README

Run these commands from the repository root, with DISCO and the example
dependencies installed.

To rerun both analyses, replace their previous results, regenerate the plots,
and render the README:

```sh
Rscript analysis/rebuild_readme.R
```

To refresh the plots and render the README using the current saved results,
without rerunning DISCO:

```sh
Rscript analysis/render_readme.R
```

Both examples write to their fixed `results/` directories. IBD input files in
`results/inputs/` are preserved. IBD figures, sample counts, pathway counts,
and significance percentages in the README come from the current saved results.
The diabetes analysis recomputes its 1,000 permutations on every run.

The user-facing README explains how to run the examples; these maintainer
commands are kept here.

After rendering, refresh the documentation website if needed:

```r
pkgdown::build_site(preview = FALSE)
```

The `Build documentation` GitHub Actions workflow also builds the site on
pull requests and pushes to `main`, and saves a `documentation-site` artifact.
Website publishing is configured separately.
