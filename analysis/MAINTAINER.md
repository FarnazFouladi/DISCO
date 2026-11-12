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

The `Build and publish documentation` GitHub Actions workflow builds the site
on pull requests and saves a downloadable `documentation-site` artifact.
Pushes to `main` also publish the site to GitHub Pages after the build succeeds.
A manual workflow run publishes only when run from `main`.

For the initial setup, select **Settings > Pages > Build and deployment >
Source > GitHub Actions** in the GitHub repository. The site address is
<https://farnazfouladi.github.io/DISCO/>.
