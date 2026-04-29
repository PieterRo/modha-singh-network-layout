# Modha-Singh Network Layout

Standalone MATLAB project for building a compact boxed 2D layout of the Modha & Singh / CoCoMac macaque long-distance connectivity network.

## Files

- `make_modha_mds_box_layout.m`: main MATLAB function
- `region_constraints.csv`: optional soft anatomical anchors for side-view orientation
- `sd01.txt`: node list
- `sd02.txt`: directed connectivity edge list
- `sd03.txt`: hierarchy / mapping edge list

Generated outputs are written to `results/` and are ignored by git by default.

## Basic Usage

```matlab
opts = struct();
opts.layout = 'grid';
opts.colorByMajorGroup = true;
opts.useRegionConstraints = true;
opts.useRefinement = true;
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

## Soft Anatomical Orientation

The project can also apply weak anterior-posterior and dorsal-ventral anchor
preferences during final box assignment without recomputing the MDS embedding.
This helps push posterior visual areas left, frontal areas right, dorsal areas
up, and ventral areas down while still letting neighborhood structure come from
connectivity.

Enable this with:

```matlab
opts.useRegionConstraints = true;
opts.constraintLambdaAP = 0.25;
opts.constraintLambdaDV = 0.20;
```

By default the function looks for `region_constraints.csv` next to
`make_modha_mds_box_layout.m`. You can point to a different file with:

```matlab
opts.regionConstraintsFile = 'my_constraints.csv';
```

The exported coordinate table includes:

- `is_anatomically_anchored`
- `anchor_ap_pref`
- `anchor_dv_pref`
- `anchor_ap_weight`
- `anchor_dv_weight`
- `anchor_source`
- `anchor_notes`

## MDS Cache

The expensive MDS embedding is cached in a MAT file so recoloring, redrawing,
or changing soft anatomical anchor weights does not need to recompute the
embedding.

Useful options:

```matlab
opts.reuseMDSCache = true;
opts.saveMDSCache = true;
opts.forceRecomputeMDS = false;
opts.mdsCacheFile = 'results/modha_boxed_layout_mds_cache.mat';
```

## Refinement Stage

The function can optionally refine the cached CoCoMac MDS coordinates before
the final box assignment. This refinement combines:

- profile preservation from the original CoCoMac MDS
- attraction between directly connected CoCoMac pairs
- optional anatomical anchor pulls
- optional Markov 2014 cortical edge weights when provided

Enable it with:

```matlab
opts.useRefinement = true;
opts.refineProfileWeight = 1.0;
opts.refineProfileK = 12;
opts.refineEdgeWeight = 0.08;
opts.refineAnatomyWeight = 0.20;
```

Useful tuning options:

```matlab
opts.refineMaxIter = 160;
opts.refineMinIter = 10;
opts.refineStepSize = 0.12;
opts.refineTolerance = 1e-4;
opts.refineEdgeTargetFraction = 0.60;
```

The exported coordinate table includes both the raw and refined coordinates:

- `mds_x`, `mds_y`
- `refined_x`, `refined_y`

## Optional Markov Weights

If you have a cortical edge-weight file derived from Markov et al. 2014, pass
it in with:

```matlab
opts.markovWeightsFile = 'markov_edge_weights.csv';
opts.markovWeightScale = 2.0;
opts.markovWeightTransform = 'auto';   % 'auto', 'log10', or 'none'
```

Expected CSV columns:

```text
source_acronym,target_acronym,weight
V1,V2,0.032
V2,V1,0.028
FEF,46d,0.011
```

The acronyms in this file should match the Modha labels used by `sd01.txt`.
The loader symmetrizes directed rows internally for the 2D refinement.

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
