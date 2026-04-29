# Modha-Singh Network Layout

Standalone MATLAB project for building a compact boxed 2D layout of the Modha & Singh / CoCoMac macaque long-distance connectivity network.

## Files

- `make_modha_mds_box_layout.m`: main MATLAB function
- `sd01.txt`: node list
- `sd02.txt`: directed connectivity edge list
- `sd03.txt`: hierarchy / mapping edge list

Generated outputs are written to `results/` and are ignored by git by default.

## Basic Usage

```matlab
opts = struct();
opts.layout = 'grid';
opts.colorByMajorGroup = true;
opts.reuseMDSCache = true;
opts.saveMDSCache = true;

OUT = make_modha_mds_box_layout( ...
    'sd01.txt', ...
    'sd02.txt', ...
    'sd03.txt', ...
    'results/modha_boxed_layout', ...
    opts);
```

## Outputs

For an output prefix such as `results/modha_boxed_layout`, the function writes:

- `results/modha_boxed_layout.png`
- `results/modha_boxed_layout.pdf`
- `results/modha_boxed_layout.svg` when SVG export succeeds
- `results/modha_boxed_layout_coordinates.csv`
- `results/modha_boxed_layout_mds_cache.mat`

## Major Groups

The function can color nodes by hierarchy-derived major group:

- `cortical`
- `thalamic`
- `basal_ganglia`
- `diencephalon`
- `lower`

Enable this with:

```matlab
opts.colorByMajorGroup = true;
```

## MDS Cache

The expensive MDS embedding is cached in a MAT file so recoloring or redrawing does not need to recompute the embedding.

Useful options:

```matlab
opts.reuseMDSCache = true;
opts.saveMDSCache = true;
opts.forceRecomputeMDS = false;
opts.mdsCacheFile = 'results/modha_boxed_layout_mds_cache.mat';
```

## Publishing To GitHub

This folder is intended to be usable as a separate git repository from the parent `Analysis_c` project.

If you already have a GitHub repo created, from this folder run:

```bash
git remote add origin <your-github-repo-url>
git branch -M main
git push -u origin main
```

If you want GitHub to create the repo for you and you have the GitHub CLI authenticated:

```bash
gh repo create <repo-name> --private --source=. --remote=origin --push
```
