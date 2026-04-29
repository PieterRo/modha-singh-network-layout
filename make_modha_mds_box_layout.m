function OUT = make_modha_mds_box_layout(sd01file, sd02file, varargin)
%MAKE_MODHA_MDS_BOX_LAYOUT Build a compact boxed MDS layout for macaque connectivity.
%
% Usage
%   OUT = make_modha_mds_box_layout(sd01file, sd02file, outPrefix)
%   OUT = make_modha_mds_box_layout(sd01file, sd02file, sd03file, outPrefix)
%   OUT = make_modha_mds_box_layout(..., opts)
%
% Required inputs
%   sd01file   Node list text file:   index acronym
%   sd02file   Edge list text file:   source_index target_index
%
% Optional inputs
%   sd03file   Hierarchy file:        parent_index child_index
%   outPrefix  Output prefix for figure/table files
%   opts       Struct with fields:
%              .useConnectedOnly      (default false)
%              .connectedMode         (default 'nonzeroDegree')
%              .layout                (default 'grid', or 'circle')
%              .aspectRatio           (default 1.4)
%              .boxWidth              (default 1)
%              .boxHeight             (default 0.45)
%              .boxGapX               (default 0.08)
%              .boxGapY               (default 0.08)
%              .fontSize              (default 6)
%              .colorByMajorGroup     (default false)
%              .useRegionConstraints  (default false)
%              .regionConstraintsFile (default local region_constraints.csv)
%              .constraintLambdaAP    (default 0.25)
%              .constraintLambdaDV    (default 0.20)
%              .useRefinement         (default false)
%              .refineProfileWeight   (default 1.0)
%              .refineProfileK        (default 12)
%              .refineEdgeWeight      (default 0.08)
%              .refineEdgeTargetFraction (default 0.60)
%              .refineAnatomyWeight   (default 0.20)
%              .refineMaxIter         (default 160)
%              .refineMinIter         (default 10)
%              .refineStepSize        (default 0.12)
%              .refineTolerance       (default 1e-4)
%              .markovWeightsFile     (default '')
%              .markovWeightScale     (default 2.0)
%              .markovWeightTransform (default 'auto')
%              .drawEdges             (default false)
%              .maxEdgesToDraw        (default 500)
%              .mdsCriterion          (default 'stress')
%              .mdsReplicates         (default 3)
%              .mdsDisplay            (default 'final')
%              .reuseMDSCache         (default true)
%              .saveMDSCache          (default true)
%              .forceRecomputeMDS     (default false)
%              .mdsCacheFile          (default [outPrefix '_mds_cache.mat'])
%              .nAngleSamples         (default 12)
%              .candidateOversubscription (default 1.10)
%              .figureVisible         (default 'off')
%              .exportSvg             (default true)
%
% Notes
%   When opts.useConnectedOnly is true, the default connected-node rule is
%   to keep nodes with nonzero in-degree or out-degree. This avoids
%   degenerate all-zero connection profiles while keeping the logic simple.
%
% Output
%   OUT is a struct with the adjacency matrix, MDS coordinates, final box
%   layout, output file paths, and the exported coordinate table.

[sd03file, outPrefix, opts] = parse_inputs(varargin{:});
opts = apply_default_opts(opts);

sd01file = to_char(sd01file);
sd02file = to_char(sd02file);
sd03file = to_char(sd03file);
outPrefix = to_char(outPrefix);

assert_file_exists(sd01file, 'sd01file');
assert_file_exists(sd02file, 'sd02file');
if ~isempty(sd03file)
    assert_file_exists(sd03file, 'sd03file');
end
ensure_output_folder(outPrefix);

[nodeIndex, labels] = read_names_list(sd01file);
nNodes = numel(nodeIndex);

connEdges = read_edge_list(sd02file, nNodes, 'connectivity');
if ~isempty(sd03file)
    mapEdges = read_edge_list(sd03file, nNodes, 'mapping');
else
    mapEdges = zeros(0, 2);
end

A = false(nNodes, nNodes);
if ~isempty(connEdges)
    A(sub2ind([nNodes, nNodes], connEdges(:,1), connEdges(:,2))) = true;
end

inDegree = full(sum(A, 1))';
outDegree = full(sum(A, 2));

if ~isempty(mapEdges)
    [parentOfAll, rootIndexAll, depthAll] = summarize_hierarchy(mapEdges, nNodes);
    rootLabelAll = repmat({''}, nNodes, 1);
    rootLabelAll(:) = labels(rootIndexAll);
    majorGroupAll = classify_major_groups(labels, parentOfAll);
else
    parentOfAll = zeros(nNodes, 1);
    rootIndexAll = nan(nNodes, 1);
    depthAll = nan(nNodes, 1);
    rootLabelAll = repmat({''}, nNodes, 1);
    majorGroupAll = repmat({'unknown'}, nNodes, 1);
end

keepMask = select_nodes(A, opts);
plotIdx = find(keepMask);

if isempty(plotIdx)
    error('make_modha_mds_box_layout:NoNodesSelected', ...
        'No nodes were selected for layout. Check opts.useConnectedOnly / opts.connectedMode.');
end

Aplot = A(plotIdx, plotIdx);
labelsPlot = labels(plotIdx);
origIndexPlot = nodeIndex(plotIdx);
majorGroupPlot = majorGroupAll(plotIdx);
constraintSpecAll = load_region_constraints(labels, opts);
constraintSpecPlot = subset_region_constraints(constraintSpecAll, plotIdx);
markovSpecPlot = load_markov_weights(labelsPlot, opts);

cacheFile = resolve_mds_cache_file(outPrefix, opts);
cacheSignature = build_mds_cache_signature(sd01file, sd02file, sd03file, plotIdx, opts);
[cacheHit, cachePayload, cacheInfo] = try_load_mds_cache(cacheFile, cacheSignature, opts);

if cacheHit
    Yplot = cachePayload.Yplot;
    stress = cachePayload.stress;
    distInfo = cachePayload.distInfo;
    mdsInfo = cachePayload.mdsInfo;
else
    X = double([Aplot Aplot']);
    [D, Dsquare, distInfo] = compute_profile_jaccard_distance(X);
    [Yplot, stress, mdsInfo] = compute_mds_layout(D, Dsquare, opts);
    cacheInfo = save_mds_cache(cacheFile, cacheSignature, Yplot, stress, distInfo, mdsInfo, opts, cacheInfo);
end
[Ylayout, refineInfo] = refine_layout_coordinates(Yplot, Aplot, labelsPlot, ...
    constraintSpecPlot, markovSpecPlot, opts);
[centerXY, gridRC, layoutInfo] = assign_box_layout(Ylayout, opts, constraintSpecPlot);
layoutInfo.refinement = refineInfo;

boxX = centerXY(:,1) - opts.boxWidth / 2;
boxY = centerXY(:,2) - opts.boxHeight / 2;

tbl = build_output_table(nodeIndex, labels, keepMask, Yplot, Ylayout, gridRC, centerXY, ...
    boxX, boxY, inDegree, outDegree, rootIndexAll, rootLabelAll, depthAll, ...
    majorGroupAll, constraintSpecAll);

fig = draw_box_layout(labelsPlot, majorGroupPlot, centerXY, Aplot, opts, layoutInfo);
files = export_outputs(fig, outPrefix, tbl, opts, cacheFile);

if ~opts.keepFigureOpen
    close(fig);
    figHandle = [];
else
    figHandle = fig;
end

if opts.verbose
    fprintf('Modha boxed layout: nodes=%d | plotted=%d | edges=%d\n', ...
        nNodes, numel(plotIdx), nnz(A));
    fprintf('MDS method: %s | stress=%g\n', mdsInfo.method, stress);
    fprintf('Saved figure/table with prefix: %s\n', outPrefix);
end

OUT = struct();
OUT.nodeIndex = nodeIndex;
OUT.labels = labels;
OUT.adjacency = A;
OUT.connectivityEdges = connEdges;
OUT.mappingEdges = mapEdges;
OUT.keepMask = keepMask;
OUT.plotIndex = plotIdx;
OUT.plotNodeIndex = origIndexPlot;
OUT.mdsXY = Yplot;
OUT.refinedXY = Ylayout;
OUT.boxCenterXY = centerXY;
OUT.gridRC = gridRC;
OUT.parentOf = parentOfAll;
OUT.majorGroup = majorGroupAll;
OUT.regionConstraints = constraintSpecAll;
OUT.markovWeights = markovSpecPlot;
OUT.inDegree = inDegree;
OUT.outDegree = outDegree;
OUT.stress = stress;
OUT.distanceInfo = distInfo;
OUT.mdsInfo = mdsInfo;
OUT.refinementInfo = refineInfo;
OUT.layoutInfo = layoutInfo;
OUT.cacheInfo = cacheInfo;
OUT.table = tbl;
OUT.files = files;
OUT.figureHandle = figHandle;
OUT.opts = opts;
end

function [sd03file, outPrefix, opts] = parse_inputs(varargin)
sd03file = '';
outPrefix = '';
opts = struct();

if isempty(varargin)
    error('make_modha_mds_box_layout:MissingInputs', ...
        'Expected at least an output prefix after sd01file and sd02file.');
end

args = varargin;
if isstruct(args{end})
    opts = args{end};
    args = args(1:end-1);
end

if numel(args) == 1
    outPrefix = args{1};
elseif numel(args) == 2
    sd03file = args{1};
    outPrefix = args{2};
else
    error('make_modha_mds_box_layout:BadSignature', ...
        ['Use make_modha_mds_box_layout(sd01, sd02, outPrefix), ' ...
         'make_modha_mds_box_layout(sd01, sd02, sd03, outPrefix), ' ...
         'or either form with opts as the last argument.']);
end
end

function opts = apply_default_opts(opts)
defaults = struct();
defaults.useConnectedOnly = false;
defaults.connectedMode = 'nonzeroDegree';
defaults.layout = 'grid';
defaults.aspectRatio = 1.4;
defaults.boxWidth = 1;
defaults.boxHeight = 0.45;
defaults.boxGapX = 0.08;
defaults.boxGapY = 0.08;
defaults.fontSize = 6;
defaults.fontName = 'Helvetica';
defaults.colorByMajorGroup = false;
defaults.useRegionConstraints = false;
defaults.regionConstraintsFile = '';
defaults.constraintLambdaAP = 0.25;
defaults.constraintLambdaDV = 0.20;
defaults.useRefinement = false;
defaults.refineProfileWeight = 1.0;
defaults.refineProfileK = 12;
defaults.refineEdgeWeight = 0.08;
defaults.refineEdgeTargetFraction = 0.60;
defaults.refineAnatomyWeight = 0.20;
defaults.refineMaxIter = 160;
defaults.refineMinIter = 10;
defaults.refineStepSize = 0.12;
defaults.refineTolerance = 1e-4;
defaults.markovWeightsFile = '';
defaults.markovWeightScale = 2.0;
defaults.markovWeightTransform = 'auto';
defaults.drawEdges = false;
defaults.maxEdgesToDraw = 500;
defaults.mdsCriterion = 'stress';
defaults.mdsReplicates = 3;
defaults.mdsDisplay = 'final';
defaults.reuseMDSCache = true;
defaults.saveMDSCache = true;
defaults.forceRecomputeMDS = false;
defaults.mdsCacheFile = '';
defaults.nAngleSamples = 12;
defaults.fitFill = 0.97;
defaults.candidateOversubscription = 1.10;
defaults.circleFillFraction = pi / 4;
defaults.useMatchpairs = true;
defaults.boxCurvature = 0.08;
defaults.lineWidth = 0.5;
defaults.outerPadding = 0.20;
defaults.boxFaceColor = [1 1 1];
defaults.boxEdgeColor = [0 0 0];
defaults.edgeColor = [0.78 0.78 0.78];
defaults.edgeLineWidth = 0.20;
defaults.figureVisible = 'off';
defaults.keepFigureOpen = false;
defaults.exportSvg = true;
defaults.pngResolution = 300;
defaults.rngSeed = [];
defaults.verbose = true;

defaultFields = fieldnames(defaults);
for i = 1:numel(defaultFields)
    name = defaultFields{i};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end

opts.layout = lower(to_char(opts.layout));
opts.connectedMode = lower(to_char(opts.connectedMode));
opts.figureVisible = lower(to_char(opts.figureVisible));
opts.mdsDisplay = lower(to_char(opts.mdsDisplay));
opts.markovWeightTransform = lower(to_char(opts.markovWeightTransform));

if ~strcmp(opts.layout, 'grid') && ~strcmp(opts.layout, 'circle')
    error('make_modha_mds_box_layout:BadLayout', ...
        'opts.layout must be ''grid'' or ''circle''.');
end

if ~strcmp(opts.mdsDisplay, 'off') && ~strcmp(opts.mdsDisplay, 'final') && ~strcmp(opts.mdsDisplay, 'iter')
    error('make_modha_mds_box_layout:BadMDSDisplay', ...
        'opts.mdsDisplay must be ''off'', ''final'', or ''iter''.');
end

if ~strcmp(opts.markovWeightTransform, 'auto') && ~strcmp(opts.markovWeightTransform, 'none') ...
        && ~strcmp(opts.markovWeightTransform, 'log10')
    error('make_modha_mds_box_layout:BadMarkovTransform', ...
        'opts.markovWeightTransform must be ''auto'', ''none'', or ''log10''.');
end
end

function value = to_char(value)
if isempty(value)
    value = '';
elseif isstring(value)
    value = char(value);
end
end

function assert_file_exists(filePath, argName)
if exist(filePath, 'file') ~= 2
    error('make_modha_mds_box_layout:MissingFile', ...
        '%s does not exist: %s', argName, filePath);
end
end

function ensure_output_folder(outPrefix)
outDir = fileparts(outPrefix);
if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
end

function cacheFile = resolve_mds_cache_file(outPrefix, opts)
if ~isempty(opts.mdsCacheFile)
    cacheFile = to_char(opts.mdsCacheFile);
else
    cacheFile = [outPrefix '_mds_cache.mat'];
end
end

function constraintSpec = load_region_constraints(labels, opts)
nNodes = numel(labels);
constraintSpec = empty_region_constraints(nNodes);
constraintSpec.file = resolve_region_constraints_file(opts);

if ~opts.useRegionConstraints
    constraintSpec.info.reason = 'disabled';
    return;
end

if isempty(constraintSpec.file) || exist(constraintSpec.file, 'file') ~= 2
    constraintSpec.info.reason = 'file_missing';
    warning('make_modha_mds_box_layout:RegionConstraintsMissing', ...
        ['opts.useRegionConstraints is true, but no region constraints file was found. ' ...
         'Continuing without anatomical anchors.']);
    return;
end

T = readtable(constraintSpec.file, 'Delimiter', ',');
requiredVars = {'acronym', 'ap_pref', 'dv_pref', 'ap_weight', 'dv_weight'};
for i = 1:numel(requiredVars)
    if ~any(strcmpi(T.Properties.VariableNames, requiredVars{i}))
        error('make_modha_mds_box_layout:BadRegionConstraintsFile', ...
            'Region constraints file is missing required column "%s": %s', ...
            requiredVars{i}, constraintSpec.file);
    end
end

constraintLabels = table_text_column_to_cellstr(get_table_var(T, 'acronym'));
apPref = double(get_table_var(T, 'ap_pref'));
dvPref = double(get_table_var(T, 'dv_pref'));
apWeight = double(get_table_var(T, 'ap_weight'));
dvWeight = double(get_table_var(T, 'dv_weight'));

if any(~isfinite(apWeight)) || any(apWeight < 0) || any(~isfinite(dvWeight)) || any(dvWeight < 0)
    error('make_modha_mds_box_layout:BadRegionConstraintsWeights', ...
        'ap_weight and dv_weight must be finite and non-negative in %s.', constraintSpec.file);
end

source = repmat({''}, height(T), 1);
notes = repmat({''}, height(T), 1);
if any(strcmpi(T.Properties.VariableNames, 'source'))
    source = table_text_column_to_cellstr(get_table_var(T, 'source'));
end
if any(strcmpi(T.Properties.VariableNames, 'notes'))
    notes = table_text_column_to_cellstr(get_table_var(T, 'notes'));
end

[tf, loc] = ismember(constraintLabels, labels);
if numel(unique(constraintLabels(tf))) ~= nnz(tf)
    error('make_modha_mds_box_layout:DuplicateRegionConstraintLabels', ...
        'Duplicate acronym entries were found among the matched rows in %s.', constraintSpec.file);
end

matchedRows = find(tf);
matchedIdx = loc(tf);
constraintSpec.apPref(matchedIdx) = apPref(matchedRows);
constraintSpec.dvPref(matchedIdx) = dvPref(matchedRows);
constraintSpec.apWeight(matchedIdx) = apWeight(matchedRows);
constraintSpec.dvWeight(matchedIdx) = dvWeight(matchedRows);
constraintSpec.source(matchedIdx) = source(matchedRows);
constraintSpec.notes(matchedIdx) = notes(matchedRows);
constraintSpec.isAnchored = ...
    (isfinite(constraintSpec.apPref) & (constraintSpec.apWeight > 0)) | ...
    (isfinite(constraintSpec.dvPref) & (constraintSpec.dvWeight > 0));
constraintSpec.active = any(constraintSpec.isAnchored);
constraintSpec.info.reason = 'loaded';
constraintSpec.info.nRows = height(T);
constraintSpec.info.nMatchedRows = numel(matchedRows);
constraintSpec.info.nUnmatchedRows = height(T) - numel(matchedRows);
constraintSpec.info.unmatchedAcronyms = constraintLabels(~tf);
constraintSpec.info.nAnchoredNodes = nnz(constraintSpec.isAnchored);

if opts.verbose
    fprintf(['Loaded region constraints from %s ' ...
        '(%d matched rows, %d anchored nodes).\n'], ...
        constraintSpec.file, constraintSpec.info.nMatchedRows, constraintSpec.info.nAnchoredNodes);
    if constraintSpec.info.nUnmatchedRows > 0
        fprintf('Region constraints ignored %d unmatched rows.\n', constraintSpec.info.nUnmatchedRows);
    end
end
end

function constraintSpec = empty_region_constraints(nNodes)
constraintSpec = struct();
constraintSpec.file = '';
constraintSpec.active = false;
constraintSpec.isAnchored = false(nNodes, 1);
constraintSpec.apPref = nan(nNodes, 1);
constraintSpec.dvPref = nan(nNodes, 1);
constraintSpec.apWeight = zeros(nNodes, 1);
constraintSpec.dvWeight = zeros(nNodes, 1);
constraintSpec.source = repmat({''}, nNodes, 1);
constraintSpec.notes = repmat({''}, nNodes, 1);
constraintSpec.info = struct('reason', '', 'nRows', 0, 'nMatchedRows', 0, ...
    'nUnmatchedRows', 0, 'unmatchedAcronyms', {{}}, 'nAnchoredNodes', 0);
end

function filePath = resolve_region_constraints_file(opts)
if ~isempty(opts.regionConstraintsFile)
    filePath = to_char(opts.regionConstraintsFile);
    return;
end

functionDir = fileparts(mfilename('fullpath'));
candidateFile = fullfile(functionDir, 'region_constraints.csv');
if exist(candidateFile, 'file') == 2
    filePath = candidateFile;
else
    filePath = '';
end
end

function value = get_table_var(T, varName)
idx = find(strcmpi(T.Properties.VariableNames, varName), 1, 'first');
if isempty(idx)
    error('make_modha_mds_box_layout:MissingTableVariable', ...
        'Table variable not found: %s', varName);
end
value = T.(T.Properties.VariableNames{idx});
end

function out = table_text_column_to_cellstr(col)
if iscellstr(col)
    out = col(:);
elseif isstring(col)
    out = cellstr(col(:));
elseif ischar(col)
    out = cellstr(col);
elseif iscell(col)
    out = cellfun(@to_char, col(:), 'UniformOutput', false);
else
    error('make_modha_mds_box_layout:BadTextColumn', ...
        'Expected a text column in region_constraints.csv.');
end
end

function subset = subset_region_constraints(constraintSpec, idx)
subset = constraintSpec;
fieldsToSubset = {'isAnchored', 'apPref', 'dvPref', 'apWeight', 'dvWeight', 'source', 'notes'};
for i = 1:numel(fieldsToSubset)
    fieldName = fieldsToSubset{i};
    subset.(fieldName) = constraintSpec.(fieldName)(idx);
end
subset.active = any(subset.isAnchored);
subset.info.nAnchoredNodes = nnz(subset.isAnchored);
end

function markovSpec = load_markov_weights(labelsPlot, opts)
nNodes = numel(labelsPlot);
markovSpec = empty_markov_weights(nNodes);
markovSpec.file = to_char(opts.markovWeightsFile);

if isempty(markovSpec.file)
    markovSpec.info.reason = 'not_provided';
    return;
end

if exist(markovSpec.file, 'file') ~= 2
    warning('make_modha_mds_box_layout:MarkovWeightsMissing', ...
        'opts.markovWeightsFile does not exist: %s. Continuing without Markov weights.', ...
        markovSpec.file);
    markovSpec.info.reason = 'file_missing';
    return;
end

T = readtable(markovSpec.file, 'Delimiter', ',');
requiredVars = {'source_acronym', 'target_acronym', 'weight'};
for i = 1:numel(requiredVars)
    if ~any(strcmpi(T.Properties.VariableNames, requiredVars{i}))
        error('make_modha_mds_box_layout:BadMarkovWeightsFile', ...
            'Markov weights file is missing required column "%s": %s', ...
            requiredVars{i}, markovSpec.file);
    end
end

srcLabel = table_text_column_to_cellstr(get_table_var(T, 'source_acronym'));
dstLabel = table_text_column_to_cellstr(get_table_var(T, 'target_acronym'));
rawWeight = double(get_table_var(T, 'weight'));
valid = isfinite(rawWeight) & (rawWeight > 0);
srcLabel = srcLabel(valid);
dstLabel = dstLabel(valid);
rawWeight = rawWeight(valid);

[srcFound, srcIdx] = ismember(srcLabel, labelsPlot);
[dstFound, dstIdx] = ismember(dstLabel, labelsPlot);
matched = srcFound & dstFound;

if ~any(matched)
    markovSpec.info.reason = 'no_matching_rows';
    warning('make_modha_mds_box_layout:MarkovWeightsNoMatches', ...
        'No rows in %s matched the plotted labels. Continuing without Markov weights.', ...
        markovSpec.file);
    return;
end

srcIdx = srcIdx(matched);
dstIdx = dstIdx(matched);
rawWeight = rawWeight(matched);
srcLabel = srcLabel(matched);
dstLabel = dstLabel(matched);

Wdir = nan(nNodes, nNodes);
for i = 1:numel(rawWeight)
    Wdir(srcIdx(i), dstIdx(i)) = rawWeight(i);
end

Wsym = nan(nNodes, nNodes);
for i = 1:nNodes
    for j = (i + 1):nNodes
        wij = Wdir(i, j);
        wji = Wdir(j, i);
        if isfinite(wij) && isfinite(wji)
            w = sqrt(wij * wji);
        elseif isfinite(wij)
            w = wij;
        elseif isfinite(wji)
            w = wji;
        else
            continue;
        end
        Wsym(i, j) = w;
        Wsym(j, i) = w;
    end
end

Wnorm = normalize_markov_weight_matrix(Wsym, opts.markovWeightTransform);
markovSpec.rawSymWeight = Wsym;
markovSpec.normSymWeight = Wnorm;
markovSpec.active = any(isfinite(Wnorm(:)) & (Wnorm(:) > 0));
markovSpec.info.reason = 'loaded';
markovSpec.info.nMatchedRows = numel(rawWeight);
markovSpec.info.nSymmetricPairs = nnz(triu(isfinite(Wnorm) & (Wnorm > 0), 1));
markovSpec.info.labelsMatched = unique([srcLabel; dstLabel]);

if opts.verbose
    fprintf(['Loaded Markov weights from %s ' ...
        '(%d matched directed rows, %d cortical pairs).\n'], ...
        markovSpec.file, markovSpec.info.nMatchedRows, markovSpec.info.nSymmetricPairs);
end
end

function markovSpec = empty_markov_weights(nNodes)
markovSpec = struct();
markovSpec.file = '';
markovSpec.active = false;
markovSpec.rawSymWeight = nan(nNodes, nNodes);
markovSpec.normSymWeight = nan(nNodes, nNodes);
markovSpec.info = struct('reason', '', 'nMatchedRows', 0, ...
    'nSymmetricPairs', 0, 'labelsMatched', {{}});
end

function Wnorm = normalize_markov_weight_matrix(Wsym, transformMode)
Wnorm = nan(size(Wsym));
positive = isfinite(Wsym) & (Wsym > 0);
if ~any(positive(:))
    return;
end

raw = Wsym(positive);
switch transformMode
    case 'none'
        values = raw;
    case 'log10'
        values = log10(raw);
    otherwise
        if max(raw) <= 1 && min(raw) < 1e-2
            values = log10(raw);
        else
            values = raw;
        end
end

vMin = min(values);
vMax = max(values);
if ~isfinite(vMin) || ~isfinite(vMax)
    return;
end
if abs(vMax - vMin) < eps
    valuesNorm = ones(size(values));
else
    valuesNorm = (values - vMin) ./ (vMax - vMin);
end

Wnorm(positive) = valuesNorm;
end

function sig = build_mds_cache_signature(sd01file, sd02file, sd03file, plotIdx, opts)
sig = struct();
sig.version = 1;
sig.distanceMethod = 'custom_binary_jaccard_v1';
sig.sd01 = file_signature(sd01file);
sig.sd02 = file_signature(sd02file);
sig.sd03 = file_signature(sd03file);
sig.plotIdx = plotIdx(:);
sig.useConnectedOnly = opts.useConnectedOnly;
sig.connectedMode = opts.connectedMode;
sig.mdsCriterion = opts.mdsCriterion;
sig.mdsReplicates = opts.mdsReplicates;
sig.rngSeed = opts.rngSeed;
end

function sig = file_signature(filePath)
if isempty(filePath)
    sig = struct('path', '', 'bytes', 0, 'datenum', NaN);
    return;
end

info = dir(filePath);
if isempty(info)
    error('make_modha_mds_box_layout:MissingFile', ...
        'Cannot build cache signature because file is missing: %s', filePath);
end

sig = struct();
sig.path = filePath;
sig.bytes = info(1).bytes;
sig.datenum = info(1).datenum;
end

function [cacheHit, payload, info] = try_load_mds_cache(cacheFile, cacheSignature, opts)
cacheHit = false;
payload = struct();
info = struct('file', cacheFile, 'used', false, 'loaded', false, 'saved', false, 'reason', '');

if opts.forceRecomputeMDS
    info.reason = 'force_recompute';
    if opts.verbose
        fprintf('Skipping cached MDS because opts.forceRecomputeMDS is true.\n');
    end
    return;
end

if ~opts.reuseMDSCache
    info.reason = 'cache_reuse_disabled';
    return;
end

if exist(cacheFile, 'file') ~= 2
    info.reason = 'cache_missing';
    return;
end

S = load(cacheFile, 'CACHE');
if ~isfield(S, 'CACHE')
    info.reason = 'cache_missing_variable';
    return;
end

CACHE = S.CACHE;
if ~isfield(CACHE, 'signature') || ~isequaln(CACHE.signature, cacheSignature)
    info.reason = 'cache_signature_mismatch';
    if opts.verbose
        fprintf('Ignoring cached MDS because the cache signature does not match current inputs.\n');
    end
    return;
end

requiredFields = {'Yplot', 'stress', 'distInfo', 'mdsInfo'};
for i = 1:numel(requiredFields)
    if ~isfield(CACHE, requiredFields{i})
        info.reason = 'cache_incomplete';
        return;
    end
end

payload = struct();
payload.Yplot = CACHE.Yplot;
payload.stress = CACHE.stress;
payload.distInfo = CACHE.distInfo;
payload.mdsInfo = CACHE.mdsInfo;

cacheHit = true;
info.used = true;
info.loaded = true;
info.reason = 'cache_loaded';

if opts.verbose
    fprintf('Loaded cached MDS from %s\n', cacheFile);
end
end

function info = save_mds_cache(cacheFile, cacheSignature, Yplot, stress, distInfo, mdsInfo, opts, info)
if ~opts.saveMDSCache
    info.reason = 'cache_not_saved';
    return;
end

CACHE = struct();
CACHE.signature = cacheSignature;
CACHE.Yplot = Yplot;
CACHE.stress = stress;
CACHE.distInfo = distInfo;
CACHE.mdsInfo = mdsInfo;

try
    save(cacheFile, 'CACHE', '-v7');
    info.saved = true;
    info.file = cacheFile;
    if isempty(info.reason) || strcmp(info.reason, 'cache_missing')
        info.reason = 'cache_saved';
    else
        info.reason = [info.reason '_and_saved'];
    end
    if opts.verbose
        fprintf('Saved MDS cache to %s\n', cacheFile);
    end
catch ME
    warning('make_modha_mds_box_layout:MDSCacheSaveFailed', ...
        'Could not save MDS cache (%s).', ME.message);
    info.saved = false;
    if isempty(info.reason)
        info.reason = 'cache_save_failed';
    end
end
end

function [nodeIndex, labels] = read_names_list(filePath)
fid = fopen(filePath, 'r');
if fid < 0
    error('make_modha_mds_box_layout:OpenFailed', ...
        'Could not open names file: %s', filePath);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
C = textscan(fid, '%f%s', ...
    'Delimiter', {' ', sprintf('\t'), ','}, ...
    'MultipleDelimsAsOne', true, ...
    'CommentStyle', '#');

idx = C{1};
rawLabels = C{2};
if isempty(idx)
    error('make_modha_mds_box_layout:EmptyNames', ...
        'No node entries were read from %s.', filePath);
end

idx = round(idx(:));
rawLabels = rawLabels(:);

if numel(unique(idx)) ~= numel(idx)
    error('make_modha_mds_box_layout:DuplicateNodeIndex', ...
        'Duplicate node indices were found in %s.', filePath);
end

nNodes = max(idx);
nodeIndex = (1:nNodes)';
labels = repmat({''}, nNodes, 1);
for i = 1:numel(idx)
    labels{idx(i)} = rawLabels{i};
end

missing = find(cellfun('isempty', labels));
if ~isempty(missing)
    error('make_modha_mds_box_layout:MissingNodeNames', ...
        'Names list is missing labels for indices: %s', num2str(missing(:)'));
end
end

function edges = read_edge_list(filePath, nNodes, listName)
fid = fopen(filePath, 'r');
if fid < 0
    error('make_modha_mds_box_layout:OpenFailed', ...
        'Could not open %s file: %s', listName, filePath);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
C = textscan(fid, '%f%f', ...
    'Delimiter', {' ', sprintf('\t'), ','}, ...
    'MultipleDelimsAsOne', true, ...
    'CommentStyle', '#');

src = C{1};
dst = C{2};
if isempty(src)
    edges = zeros(0, 2);
    return;
end

edges = [round(src(:)), round(dst(:))];
bad = any(~isfinite(edges), 2);
edges = edges(~bad, :);

if any(edges(:) < 1) || any(edges(:) > nNodes)
    error('make_modha_mds_box_layout:EdgeIndexOutOfRange', ...
        'At least one %s edge references a node outside 1..%d.', listName, nNodes);
end

edges = unique(edges, 'rows');
end

function [parentOf, rootIndex, depth] = summarize_hierarchy(mapEdges, nNodes)
parentOf = zeros(nNodes, 1);

for i = 1:size(mapEdges, 1)
    p = mapEdges(i, 1);
    c = mapEdges(i, 2);
    if parentOf(c) ~= 0 && parentOf(c) ~= p
        error('make_modha_mds_box_layout:MultipleParents', ...
            'Node %d has multiple parents in sd03.', c);
    end
    parentOf(c) = p;
end

rootIndex = zeros(nNodes, 1);
depth = zeros(nNodes, 1);

for i = 1:nNodes
    cur = i;
    d = 0;
    seen = false(nNodes, 1);
    while parentOf(cur) ~= 0
        if seen(cur)
            error('make_modha_mds_box_layout:HierarchyCycle', ...
                'Cycle detected in sd03 hierarchy near node %d.', cur);
        end
        seen(cur) = true;
        cur = parentOf(cur);
        d = d + 1;
    end
    rootIndex(i) = cur;
    depth(i) = d;
end
end

function majorGroup = classify_major_groups(labels, parentOf)
nNodes = numel(labels);
majorGroup = repmat({'lower'}, nNodes, 1);

if isempty(parentOf)
    majorGroup(:) = {'unknown'};
    return;
end

for i = 1:nNodes
    tokens = ancestor_tokens(i, labels, parentOf);

    if has_ancestor_match(tokens, {'tha', 'thalamus'}, {'thalam'})
        majorGroup{i} = 'thalamic';
    elseif has_ancestor_match(tokens, ...
            {'bg', 'basalganglia', 'basalgangliaaccordingtogmdefinition'}, ...
            {'basalganglia'})
        majorGroup{i} = 'basal_ganglia';
    elseif has_ancestor_match(tokens, ...
            {'gmcerebralcortex', 'cerebralcortex', 'cortex', 'cx'}, ...
            {'cortex', 'cortical', 'frontallobe', 'parietallobe', ...
             'temporallobe', 'occipitallobe', 'cingulate', 'insula'})
        majorGroup{i} = 'cortical';
    elseif has_ancestor_match(tokens, {'die', 'diencephalon', 'hyp', 'hypothalamus'}, ...
            {'dienceph', 'hypothalam'})
        majorGroup{i} = 'diencephalon';
    else
        majorGroup{i} = 'lower';
    end
end
end

function tokens = ancestor_tokens(nodeIdx, labels, parentOf)
tokens = cell(0, 1);
cur = nodeIdx;
seen = false(numel(labels), 1);

while cur > 0
    if seen(cur)
        break;
    end
    seen(cur) = true;
    tokens{end+1,1} = normalize_token(labels{cur}); %#ok<AGROW>
    cur = parentOf(cur);
end
end

function tf = has_ancestor_match(tokens, exactTokens, containsTokens)
tf = false;

for i = 1:numel(tokens)
    tok = tokens{i};
    if any(strcmp(tok, exactTokens))
        tf = true;
        return;
    end

    for j = 1:numel(containsTokens)
        if contains(tok, containsTokens{j})
            tf = true;
            return;
        end
    end
end
end

function tok = normalize_token(str)
tok = lower(to_char(str));
tok = regexprep(tok, '[^a-z0-9]', '');
end

function keepMask = select_nodes(A, opts)
nNodes = size(A, 1);
keepMask = true(nNodes, 1);

if ~opts.useConnectedOnly
    return;
end

deg = full(sum(A, 1))' + full(sum(A, 2));
switch opts.connectedMode
    case 'nonzerodegree'
        keepMask = deg > 0;

    case 'largestweakcomponent'
        if nnz(A) == 0
            keepMask = false(nNodes, 1);
            return;
        end
        G = graph(double(A | A'));
        bins = conncomp(G);
        counts = accumarray(bins(:), 1);
        [~, bestBin] = max(counts);
        keepMask = (bins(:) == bestBin) & (deg > 0);

    otherwise
        error('make_modha_mds_box_layout:BadConnectedMode', ...
            'Unknown opts.connectedMode: %s', opts.connectedMode);
end
end

function [D, Dsquare, info] = compute_profile_jaccard_distance(X)
% Define Jaccard distance for all-zero profiles as 0 rather than NaN.
Xbin = double(X ~= 0);
nNodes = size(Xbin, 1);
rowSum = sum(Xbin, 2);
intersection = Xbin * Xbin';
unionCount = rowSum + rowSum' - intersection;

Dsquare = zeros(nNodes, nNodes);
hasUnion = unionCount > 0;
Dsquare(hasUnion) = 1 - intersection(hasUnion) ./ unionCount(hasUnion);
Dsquare(~hasUnion) = 0;

Dsquare(1:nNodes+1:end) = 0;
Dsquare = (Dsquare + Dsquare') / 2;

info = struct();
info.method = 'custom_binary_jaccard';
info.nZeroProfileNodes = nnz(rowSum == 0);
info.nZeroUnionPairs = nnz(triu(~hasUnion, 1));

D = squareform(Dsquare);
end

function [Y, stress, info] = compute_mds_layout(D, Dsquare, opts)
info = struct();
info.method = '';
info.criterion = opts.mdsCriterion;
info.replicates = opts.mdsReplicates;
info.display = opts.mdsDisplay;

if isempty(D) || all(D == 0)
    nNodes = size(Dsquare, 1);
    Y = zeros(nNodes, 2);
    stress = 0;
    info.method = 'degenerate_zero_distance';
    info.totalElapsedSeconds = 0;
    return;
end

if ~isempty(opts.rngSeed)
    rng(opts.rngSeed);
end

mdsOptions = statset('Display', opts.mdsDisplay);
nNodes = size(Dsquare, 1);

if opts.verbose
    fprintf(['Starting MDS: n=%d nodes | criterion=%s | replicates=%d | ' ...
        'display=%s\n'], nNodes, opts.mdsCriterion, opts.mdsReplicates, opts.mdsDisplay);
    fprintf('MDS attempt 1/2: default initialization...\n');
end

tAttempt1 = tic;
try
    [Y, stress] = mdscale(D, 2, ...
        'Criterion', opts.mdsCriterion, ...
        'Replicates', opts.mdsReplicates, ...
        'Options', mdsOptions);
    info.method = 'mdscale';
    info.attempt1Seconds = toc(tAttempt1);
    info.totalElapsedSeconds = info.attempt1Seconds;
    if opts.verbose
        fprintf('MDS attempt 1/2 finished in %.1f s.\n', info.attempt1Seconds);
    end
catch ME1
    attempt1Seconds = toc(tAttempt1);
    info.attempt1Seconds = attempt1Seconds;
    if opts.verbose
        fprintf('MDS attempt 1/2 failed after %.1f s.\n', attempt1Seconds);
    end
    warning('make_modha_mds_box_layout:MDSFallback', ...
        ['mdscale failed with default initialization (%s). ' ...
         'Retrying with random initialization.'], ME1.message);
    if opts.verbose
        fprintf('MDS attempt 2/2: random initialization...\n');
    end
    tAttempt2 = tic;
    try
        [Y, stress] = mdscale(D, 2, ...
            'Criterion', opts.mdsCriterion, ...
            'Replicates', opts.mdsReplicates, ...
            'Start', 'random', ...
            'Options', mdsOptions);
        info.method = 'mdscale_random_start';
        info.attempt2Seconds = toc(tAttempt2);
        info.totalElapsedSeconds = attempt1Seconds + info.attempt2Seconds;
        if opts.verbose
            fprintf(['MDS attempt 2/2 finished in %.1f s ' ...
                '(total %.1f s).\n'], info.attempt2Seconds, info.totalElapsedSeconds);
        end
    catch ME2
        attempt2Seconds = toc(tAttempt2);
        info.attempt2Seconds = attempt2Seconds;
        if opts.verbose
            fprintf(['MDS attempt 2/2 failed after %.1f s ' ...
                '(total %.1f s).\n'], attempt2Seconds, attempt1Seconds + attempt2Seconds);
            fprintf('Falling back to cmdscale...\n');
        end
        warning('make_modha_mds_box_layout:MDSFallback', ...
            'Random-start mdscale failed (%s). Falling back to cmdscale.', ME2.message);
        tCmd = tic;
        [Y, eigvals] = cmdscale(Dsquare, 2);
        if size(Y, 2) < 2
            Y(:,2) = 0;
        end
        if isempty(Y)
            Y = zeros(size(Dsquare, 1), 2);
        end
        stress = NaN;
        info.method = 'cmdscale';
        info.eigenvalues = eigvals;
        info.cmdscaleSeconds = toc(tCmd);
        info.totalElapsedSeconds = attempt1Seconds + attempt2Seconds + info.cmdscaleSeconds;
        if opts.verbose
            fprintf('cmdscale fallback finished in %.1f s (total %.1f s).\n', ...
                info.cmdscaleSeconds, info.totalElapsedSeconds);
        end
    end
end

if size(Y, 2) < 2
    Y(:,2) = 0;
end
end

function [Yrefined, info] = refine_layout_coordinates(Y0, Aplot, labelsPlot, ...
    constraintSpec, markovSpec, opts)
Yrefined = Y0;
info = struct();
info.enabled = opts.useRefinement;
info.reason = '';
info.iterations = 0;
info.profilePairs = 0;
info.edgePairs = 0;
info.markovPairs = 0;
info.anatomyAnchors = 0;
info.profileLoss = 0;
info.edgeLoss = 0;
info.anatomyLoss = 0;
info.totalLoss = 0;
info.delta = 0;

if ~opts.useRefinement
    info.reason = 'disabled';
    return;
end

if size(Y0, 1) < 2
    info.reason = 'too_few_nodes';
    return;
end

profileSpec = build_profile_refinement_spec(Y0, opts);
edgeSpec = build_edge_refinement_spec(Aplot, markovSpec, Y0, opts);
anatomySpec = build_anatomy_refinement_spec(constraintSpec, Y0);

info.profilePairs = numel(profileSpec.i);
info.edgePairs = numel(edgeSpec.i);
info.markovPairs = edgeSpec.nMarkovPairs;
info.anatomyAnchors = nnz(anatomySpec.isActive);

if info.profilePairs == 0 && info.edgePairs == 0 && info.anatomyAnchors == 0
    info.reason = 'no_active_terms';
    return;
end

Y = center_and_match_reference_scale(Y0, Y0);
lastLoss = inf;

for iter = 1:opts.refineMaxIter
    [forceProfile, lossProfile] = compute_profile_force(Y, profileSpec, opts.refineProfileWeight);
    [forceEdge, lossEdge] = compute_edge_force(Y, edgeSpec, opts.refineEdgeWeight);
    [forceAnatomy, lossAnatomy] = compute_anatomy_force(Y, anatomySpec, opts.refineAnatomyWeight);

    totalForce = forceProfile + forceEdge + forceAnatomy;
    Ynext = Y + opts.refineStepSize * totalForce;
    Ynext = center_and_match_reference_scale(Ynext, Y0);

    delta = sqrt(mean(sum((Ynext - Y).^2, 2)));
    totalLoss = lossProfile + lossEdge + lossAnatomy;

    Y = Ynext;
    info.iterations = iter;
    info.profileLoss = lossProfile;
    info.edgeLoss = lossEdge;
    info.anatomyLoss = lossAnatomy;
    info.totalLoss = totalLoss;
    info.delta = delta;

    if opts.verbose && (iter == 1 || mod(iter, 25) == 0 || delta < opts.refineTolerance)
        fprintf(['Refinement iter %d: total=%.5f | profile=%.5f | edge=%.5f | ' ...
            'anat=%.5f | delta=%.6f\n'], ...
            iter, totalLoss, lossProfile, lossEdge, lossAnatomy, delta);
    end

    if iter >= opts.refineMinIter && delta < opts.refineTolerance
        info.reason = 'converged';
        break;
    end
    if iter >= opts.refineMinIter && abs(lastLoss - totalLoss) < opts.refineTolerance * max(1, lastLoss)
        info.reason = 'loss_plateau';
        break;
    end
    lastLoss = totalLoss;
end

if isempty(info.reason)
    if info.iterations >= opts.refineMaxIter
        info.reason = 'max_iter';
    else
        info.reason = 'completed';
    end
end

Yrefined = Y;
info.labels = labelsPlot;
end

function spec = build_profile_refinement_spec(Y0, opts)
nNodes = size(Y0, 1);
K = min(max(1, round(opts.refineProfileK)), max(1, nNodes - 1));
spec = struct('i', zeros(0,1), 'j', zeros(0,1), 'targetDist', zeros(0,1));

if nNodes < 2 || K < 1 || opts.refineProfileWeight <= 0
    return;
end

D2 = pdist2(Y0, Y0, 'squaredeuclidean');
D2(1:nNodes+1:end) = inf;
pairList = zeros(nNodes * K, 2);
nPairs = 0;

for i = 1:nNodes
    [~, order] = sort(D2(i, :), 'ascend');
    nnIdx = order(1:K);
    pairs = sort([repmat(i, K, 1), nnIdx(:)], 2);
    pairList(nPairs + (1:K), :) = pairs;
    nPairs = nPairs + K;
end

pairList = unique(pairList(1:nPairs, :), 'rows');
keep = pairList(:,1) ~= pairList(:,2);
pairList = pairList(keep, :);

spec.i = pairList(:,1);
spec.j = pairList(:,2);
spec.targetDist = sqrt(sum((Y0(spec.i,:) - Y0(spec.j,:)).^2, 2));
end

function spec = build_edge_refinement_spec(Aplot, markovSpec, Y0, opts)
nNodes = size(Aplot, 1);
spec = struct('i', zeros(0,1), 'j', zeros(0,1), 'springWeight', zeros(0,1), ...
    'targetDist', zeros(0,1), 'nMarkovPairs', 0);

if nNodes < 2 || opts.refineEdgeWeight <= 0
    return;
end

adjUndir = Aplot | Aplot';
markovPairMask = false(nNodes, nNodes);
markovNorm = nan(nNodes, nNodes);
if markovSpec.active
    markovNorm = markovSpec.normSymWeight;
    markovPairMask = isfinite(markovNorm) & (markovNorm > 0);
end

pairMask = triu(adjUndir | markovPairMask, 1);
[ii, jj] = find(pairMask);
if isempty(ii)
    return;
end

springWeight = double(adjUndir(sub2ind([nNodes, nNodes], ii, jj)));
hasMarkov = markovPairMask(sub2ind([nNodes, nNodes], ii, jj));
springWeight(hasMarkov) = springWeight(hasMarkov) + ...
    opts.markovWeightScale * markovNorm(sub2ind([nNodes, nNodes], ii(hasMarkov), jj(hasMarkov)));

positiveWeight = springWeight > 0;
ii = ii(positiveWeight);
jj = jj(positiveWeight);
springWeight = springWeight(positiveWeight);
if isempty(ii)
    return;
end

spec.i = ii;
spec.j = jj;
spec.springWeight = springWeight;
profileSpec = build_profile_refinement_spec(Y0, opts);
if isempty(profileSpec.targetDist)
    targetScale = median(pdist(Y0));
else
    targetScale = median(profileSpec.targetDist);
end
if ~isfinite(targetScale) || targetScale <= 0
    targetScale = 1;
end
spec.targetDist = opts.refineEdgeTargetFraction * targetScale * ones(size(ii));
spec.nMarkovPairs = nnz(hasMarkov(positiveWeight));
end

function spec = build_anatomy_refinement_spec(constraintSpec, Y0)
nNodes = size(Y0, 1);
spec = struct();
spec.isActive = false(nNodes, 1);
spec.targetXY = zeros(nNodes, 2);
spec.weightXY = zeros(nNodes, 2);

if ~constraintSpec.active
    return;
end

Ycentered = Y0 - mean(Y0, 1);
scaleX = max(abs(Ycentered(:,1)));
scaleY = max(abs(Ycentered(:,2)));
if ~isfinite(scaleX) || scaleX <= 0
    scaleX = 1;
end
if ~isfinite(scaleY) || scaleY <= 0
    scaleY = 1;
end

hasAP = isfinite(constraintSpec.apPref) & (constraintSpec.apWeight > 0);
hasDV = isfinite(constraintSpec.dvPref) & (constraintSpec.dvWeight > 0);
spec.isActive = hasAP | hasDV;
spec.targetXY(hasAP, 1) = constraintSpec.apPref(hasAP) * scaleX;
spec.targetXY(hasDV, 2) = constraintSpec.dvPref(hasDV) * scaleY;
spec.weightXY(:,1) = constraintSpec.apWeight;
spec.weightXY(:,2) = constraintSpec.dvWeight;
end

function [force, loss] = compute_profile_force(Y, spec, weightScale)
force = zeros(size(Y));
loss = 0;
pairCount = zeros(size(Y,1), 1);

if isempty(spec.i) || weightScale <= 0
    return;
end

nPairs = numel(spec.i);
for p = 1:nPairs
    i = spec.i(p);
    j = spec.j(p);
    delta = Y(j,:) - Y(i,:);
    dist = hypot(delta(1), delta(2));
    if dist < eps
        dist = eps;
    end
    err = dist - spec.targetDist(p);
    unitVec = delta / dist;
    pull = weightScale * err * unitVec;
    force(i,:) = force(i,:) + pull;
    force(j,:) = force(j,:) - pull;
    pairCount(i) = pairCount(i) + 1;
    pairCount(j) = pairCount(j) + 1;
    loss = loss + weightScale * (err ^ 2) / nPairs;
end

nonzero = pairCount > 0;
force(nonzero,:) = force(nonzero,:) ./ pairCount(nonzero);
end

function [force, loss] = compute_edge_force(Y, spec, weightScale)
force = zeros(size(Y));
loss = 0;
pairCount = zeros(size(Y,1), 1);

if isempty(spec.i) || weightScale <= 0
    return;
end

nPairs = numel(spec.i);
for p = 1:nPairs
    i = spec.i(p);
    j = spec.j(p);
    delta = Y(j,:) - Y(i,:);
    dist = hypot(delta(1), delta(2));
    if dist < eps
        dist = eps;
    end
    err = dist - spec.targetDist(p);
    unitVec = delta / dist;
    localWeight = weightScale * spec.springWeight(p);
    pull = localWeight * err * unitVec;
    force(i,:) = force(i,:) + pull;
    force(j,:) = force(j,:) - pull;
    pairCount(i) = pairCount(i) + 1;
    pairCount(j) = pairCount(j) + 1;
    loss = loss + localWeight * (err ^ 2) / nPairs;
end

nonzero = pairCount > 0;
force(nonzero,:) = force(nonzero,:) ./ pairCount(nonzero);
end

function [force, loss] = compute_anatomy_force(Y, spec, weightScale)
force = zeros(size(Y));
loss = 0;

if ~any(spec.isActive) || weightScale <= 0
    return;
end

idx = find(spec.isActive);
nActive = numel(idx);
for k = 1:nActive
    i = idx(k);
    dx = spec.targetXY(i,1) - Y(i,1);
    dy = spec.targetXY(i,2) - Y(i,2);
    wx = weightScale * spec.weightXY(i,1);
    wy = weightScale * spec.weightXY(i,2);
    force(i,1) = force(i,1) + wx * dx;
    force(i,2) = force(i,2) + wy * dy;
    loss = loss + (wx * dx ^ 2 + wy * dy ^ 2) / nActive;
end
end

function Ynorm = center_and_match_reference_scale(Y, Yref)
Ynorm = Y - mean(Y, 1);
refCentered = Yref - mean(Yref, 1);
refRms = sqrt(mean(sum(refCentered .^ 2, 2)));
curRms = sqrt(mean(sum(Ynorm .^ 2, 2)));

if ~isfinite(refRms) || refRms <= 0
    refRms = 1;
end
if ~isfinite(curRms) || curRms <= 0
    curRms = 1;
end

Ynorm = Ynorm * (refRms / curRms);
end

function [centerXY, gridRC, layoutInfo] = assign_box_layout(Y, opts, constraintSpec)
nNodes = size(Y, 1);
cellW = opts.boxWidth + opts.boxGapX;
cellH = opts.boxHeight + opts.boxGapY;

shapeList = candidate_shapes(nNodes, opts, cellW, cellH);
angles = linspace(0, pi, max(1, opts.nAngleSamples) + 1);
angles(end) = [];
if isempty(angles)
    angles = 0;
end
transformList = orientation_transform_list(constraintSpec);

bestCost = inf;
bestCenterXY = [];
bestGridRC = [];
bestInfo = struct();

for s = 1:size(shapeList, 1)
    nRows = shapeList(s, 1);
    nCols = shapeList(s, 2);
    [candidateXY, candidateRC] = make_candidate_cells(nRows, nCols, nNodes, opts, cellW, cellH);
    candidateNormXY = normalize_candidate_coordinates(candidateXY);

    for t = 1:size(transformList, 1)
        reflectX = transformList(t, 1);
        reflectY = transformList(t, 2);

        for a = 1:numel(angles)
            Yrot = orient_points(Y, angles(a), reflectX, reflectY);
            Yfit = fit_points_to_candidates(Yrot, candidateXY, opts.fitFill);
            [assignIdx, totalCost, method, costInfo] = assign_points_to_cells( ...
                Yfit, candidateXY, candidateNormXY, constraintSpec, opts);

            if totalCost < bestCost
                bestCost = totalCost;
                bestCenterXY = candidateXY(assignIdx, :);
                bestGridRC = candidateRC(assignIdx, :);
                bestInfo.nRows = nRows;
                bestInfo.nCols = nCols;
                bestInfo.angleDeg = angles(a) * 180 / pi;
                bestInfo.reflectX = reflectX;
                bestInfo.reflectY = reflectY;
                bestInfo.assignmentMethod = method;
                bestInfo.totalCost = totalCost;
                bestInfo.baseCost = costInfo.baseTotal;
                bestInfo.constraintCost = costInfo.constraintTotal;
                bestInfo.nAnchoredNodes = nnz(constraintSpec.isAnchored);
            end
        end
    end
end

centerXY = bestCenterXY;
gridRC = bestGridRC;

xShift = -min(centerXY(:,1)) + opts.boxWidth / 2 + opts.outerPadding;
yShift = -min(centerXY(:,2)) + opts.boxHeight / 2 + opts.outerPadding;
centerXY(:,1) = centerXY(:,1) + xShift;
centerXY(:,2) = centerXY(:,2) + yShift;

layoutInfo = bestInfo;
layoutInfo.totalWidth = max(centerXY(:,1)) - min(centerXY(:,1)) + opts.boxWidth;
layoutInfo.totalHeight = max(centerXY(:,2)) - min(centerXY(:,2)) + opts.boxHeight;
layoutInfo.cellWidth = cellW;
layoutInfo.cellHeight = cellH;
end

function transformList = orientation_transform_list(constraintSpec)
if constraintSpec.active
    transformList = [1 1; -1 1; 1 -1; -1 -1];
else
    transformList = [1 1];
end
end

function shapeList = candidate_shapes(nNodes, opts, cellW, cellH)
switch opts.layout
    case 'grid'
        targetCells = nNodes;
    case 'circle'
        targetCells = ceil(nNodes / max(opts.circleFillFraction, 0.1));
    otherwise
        error('make_modha_mds_box_layout:BadLayout', ...
            'Unknown layout: %s', opts.layout);
end

idealCols = sqrt(targetCells * opts.aspectRatio * cellH / max(cellW, eps));
if ~isfinite(idealCols) || idealCols < 1
    idealCols = sqrt(targetCells);
end

colCandidates = unique(max(1, round(idealCols + (-3:3))));
shapeList = zeros(numel(colCandidates), 2);
nKeep = 0;

for i = 1:numel(colCandidates)
    nCols = colCandidates(i);
    nRows = ceil(targetCells / nCols);

    if strcmp(opts.layout, 'circle')
        while true
            [candidateXY, ~] = make_candidate_cells(nRows, nCols, nNodes, opts, cellW, cellH);
            if size(candidateXY, 1) >= nNodes
                break;
            end
            nRows = nRows + 1;
        end
    end

    nKeep = nKeep + 1;
    shapeList(nKeep, :) = [nRows, nCols];
end

shapeList = unique(shapeList(1:nKeep, :), 'rows');
end

function [candidateXY, candidateRC] = make_candidate_cells(nRows, nCols, nNodes, opts, cellW, cellH)
x = ((1:nCols) - (nCols + 1) / 2) * cellW;
y = ((nRows:-1:1) - (nRows + 1) / 2) * cellH;
[Xg, Yg] = meshgrid(x, y);
candidateXYAll = [Xg(:), Yg(:)];

[rowGrid, colGrid] = ndgrid(1:nRows, 1:nCols);
candidateRCAll = [rowGrid(:), colGrid(:)];

if strcmp(opts.layout, 'grid')
    keep = true(size(candidateXYAll, 1), 1);
else
    rx = max(abs(x));
    ry = max(abs(y));
    if rx == 0
        rx = cellW / 2;
    end
    if ry == 0
        ry = cellH / 2;
    end

    radial = (candidateXYAll(:,1) ./ rx).^2 + (candidateXYAll(:,2) ./ ry).^2;
    [~, order] = sort(radial, 'ascend');
    nKeep = min(numel(order), max(nNodes, ceil(nNodes * opts.candidateOversubscription)));
    keep = false(numel(order), 1);
    keep(order(1:nKeep)) = true;
end

candidateXY = candidateXYAll(keep, :);
candidateRC = candidateRCAll(keep, :);
end

function Yrot = rotate_points(Y, angleRad)
R = [cos(angleRad), -sin(angleRad); sin(angleRad), cos(angleRad)];
Yrot = Y * R';
end

function Yout = orient_points(Y, angleRad, reflectX, reflectY)
Ytmp = Y;
Ytmp(:,1) = reflectX * Ytmp(:,1);
Ytmp(:,2) = reflectY * Ytmp(:,2);
Yout = rotate_points(Ytmp, angleRad);
end

function Yfit = fit_points_to_candidates(Y, candidateXY, fitFill)
Yc = Y - mean(Y, 1);
targetCenter = mean(candidateXY, 1);

xRange = max(Yc(:,1)) - min(Yc(:,1));
yRange = max(Yc(:,2)) - min(Yc(:,2));
targetW = max(candidateXY(:,1)) - min(candidateXY(:,1));
targetH = max(candidateXY(:,2)) - min(candidateXY(:,2));

if xRange <= 0
    xRange = 1;
end
if yRange <= 0
    yRange = 1;
end
if targetW <= 0
    targetW = 1;
end
if targetH <= 0
    targetH = 1;
end

scale = fitFill * min(targetW / xRange, targetH / yRange);
if ~isfinite(scale) || scale <= 0
    scale = 1;
end

Yfit = Yc * scale;
Yfit(:,1) = Yfit(:,1) + targetCenter(1);
Yfit(:,2) = Yfit(:,2) + targetCenter(2);
end

function candidateNormXY = normalize_candidate_coordinates(candidateXY)
candidateNormXY = zeros(size(candidateXY));
candidateNormXY(:,1) = normalize_signed_coordinate(candidateXY(:,1));
candidateNormXY(:,2) = normalize_signed_coordinate(candidateXY(:,2));
end

function coordNorm = normalize_signed_coordinate(coord)
maxAbs = max(abs(coord));
if ~isfinite(maxAbs) || maxAbs <= 0
    coordNorm = zeros(size(coord));
else
    coordNorm = coord / maxAbs;
end
end

function penaltyMatrix = build_constraint_penalty(candidateNormXY, constraintSpec, opts)
nNodes = numel(constraintSpec.isAnchored);
nCells = size(candidateNormXY, 1);
penaltyMatrix = zeros(nNodes, nCells);

if ~constraintSpec.active
    return;
end

hasAP = isfinite(constraintSpec.apPref) & (constraintSpec.apWeight > 0);
if any(hasAP)
    apDelta = bsxfun(@minus, candidateNormXY(:,1)', constraintSpec.apPref(hasAP));
    penaltyMatrix(hasAP, :) = penaltyMatrix(hasAP, :) + bsxfun(@times, ...
        apDelta .^ 2, opts.constraintLambdaAP * constraintSpec.apWeight(hasAP));
end

hasDV = isfinite(constraintSpec.dvPref) & (constraintSpec.dvWeight > 0);
if any(hasDV)
    dvDelta = bsxfun(@minus, candidateNormXY(:,2)', constraintSpec.dvPref(hasDV));
    penaltyMatrix(hasDV, :) = penaltyMatrix(hasDV, :) + bsxfun(@times, ...
        dvDelta .^ 2, opts.constraintLambdaDV * constraintSpec.dvWeight(hasDV));
end
end

function [assignIdx, totalCost, method, costInfo] = assign_points_to_cells( ...
    Yfit, candidateXY, candidateNormXY, constraintSpec, opts)
baseCostMatrix = pdist2(Yfit, candidateXY, 'squaredeuclidean');
constraintCostMatrix = build_constraint_penalty(candidateNormXY, constraintSpec, opts);
costMatrix = baseCostMatrix + constraintCostMatrix;

if opts.useMatchpairs && exist('matchpairs', 'file') == 2
    unmatchedCost = max(costMatrix(:));
    if ~isfinite(unmatchedCost) || unmatchedCost <= 0
        unmatchedCost = 1;
    end
    unmatchedCost = unmatchedCost * 1e6 + 1;

    pairs = matchpairs(costMatrix, unmatchedCost);
    if size(pairs, 1) == size(costMatrix, 1)
        [~, order] = sort(pairs(:,1));
        pairs = pairs(order, :);
        assignIdx = pairs(:,2);
        linearIdx = sub2ind(size(costMatrix), (1:size(costMatrix,1))', assignIdx);
        totalCost = sum(costMatrix(linearIdx));
        costInfo = struct();
        costInfo.baseTotal = sum(baseCostMatrix(linearIdx));
        costInfo.constraintTotal = sum(constraintCostMatrix(linearIdx));
        costInfo.total = totalCost;
        method = 'matchpairs';
        return;
    end
end

[assignIdx, totalCost] = greedy_assignment(costMatrix, Yfit);
linearIdx = sub2ind(size(costMatrix), (1:size(costMatrix,1))', assignIdx);
costInfo = struct();
costInfo.baseTotal = sum(baseCostMatrix(linearIdx));
costInfo.constraintTotal = sum(constraintCostMatrix(linearIdx));
costInfo.total = totalCost;
method = 'greedy';
end

function [assignIdx, totalCost] = greedy_assignment(costMatrix, Yfit)
nNodes = size(costMatrix, 1);
nCells = size(costMatrix, 2);
available = true(1, nCells);
assignIdx = zeros(nNodes, 1);

radius2 = sum(Yfit.^2, 2);
[~, order] = sort(radius2, 'descend');

for k = 1:nNodes
    i = order(k);
    rowCost = costMatrix(i, :);
    rowCost(~available) = inf;
    [~, j] = min(rowCost);
    assignIdx(i) = j;
    available(j) = false;
end

[assignIdx, totalCost] = improve_by_swaps(assignIdx, costMatrix);
end

function [assignIdx, totalCost] = improve_by_swaps(assignIdx, costMatrix)
nNodes = numel(assignIdx);

for iter = 1:3
    improved = false;
    for i = 1:(nNodes - 1)
        ji = assignIdx(i);
        costI = costMatrix(i, ji);
        for k = (i + 1):nNodes
            jk = assignIdx(k);
            delta = costMatrix(i, jk) + costMatrix(k, ji) ...
                - (costI + costMatrix(k, jk));
            if delta < -1e-12
                assignIdx(i) = jk;
                assignIdx(k) = ji;
                ji = assignIdx(i);
                costI = costMatrix(i, ji);
                improved = true;
            end
        end
    end
    if ~improved
        break;
    end
end

totalCost = sum(costMatrix(sub2ind(size(costMatrix), (1:nNodes)', assignIdx)));
end

function tbl = build_output_table(nodeIndex, labels, keepMask, Yplot, Yrefined, gridRC, centerXY, ...
    boxX, boxY, inDegree, outDegree, rootIndexAll, rootLabelAll, depthAll, ...
    majorGroupAll, constraintSpecAll)
nNodes = numel(nodeIndex);

mdsX = nan(nNodes, 1);
mdsY = nan(nNodes, 1);
refinedX = nan(nNodes, 1);
refinedY = nan(nNodes, 1);
gridRow = nan(nNodes, 1);
gridCol = nan(nNodes, 1);
boxCenterX = nan(nNodes, 1);
boxCenterY = nan(nNodes, 1);
boxLeftX = nan(nNodes, 1);
boxBottomY = nan(nNodes, 1);
isAnchored = constraintSpecAll.isAnchored(:);
anchorApPref = constraintSpecAll.apPref(:);
anchorDvPref = constraintSpecAll.dvPref(:);
anchorApWeight = constraintSpecAll.apWeight(:);
anchorDvWeight = constraintSpecAll.dvWeight(:);
anchorSource = constraintSpecAll.source(:);
anchorNotes = constraintSpecAll.notes(:);

plotIdx = find(keepMask);
mdsX(plotIdx) = Yplot(:,1);
mdsY(plotIdx) = Yplot(:,2);
refinedX(plotIdx) = Yrefined(:,1);
refinedY(plotIdx) = Yrefined(:,2);
gridRow(plotIdx) = gridRC(:,1);
gridCol(plotIdx) = gridRC(:,2);
boxCenterX(plotIdx) = centerXY(:,1);
boxCenterY(plotIdx) = centerXY(:,2);
boxLeftX(plotIdx) = boxX;
boxBottomY(plotIdx) = boxY;

tbl = table( ...
    nodeIndex(:), ...
    labels(:), ...
    keepMask(:), ...
    mdsX(:), ...
    mdsY(:), ...
    refinedX(:), ...
    refinedY(:), ...
    gridRow(:), ...
    gridCol(:), ...
    boxCenterX(:), ...
    boxCenterY(:), ...
    boxLeftX(:), ...
    boxBottomY(:), ...
    inDegree(:), ...
    outDegree(:), ...
    rootIndexAll(:), ...
    rootLabelAll(:), ...
    majorGroupAll(:), ...
    isAnchored(:), ...
    anchorApPref(:), ...
    anchorDvPref(:), ...
    anchorApWeight(:), ...
    anchorDvWeight(:), ...
    anchorSource(:), ...
    anchorNotes(:), ...
    depthAll(:), ...
    'VariableNames', { ...
    'index', 'acronym', 'is_plotted', 'mds_x', 'mds_y', 'refined_x', 'refined_y', ...
    'grid_row', 'grid_col', ...
    'box_center_x', 'box_center_y', 'box_x', 'box_y', ...
    'in_degree', 'out_degree', 'root_index', 'root_acronym', 'major_group', ...
    'is_anatomically_anchored', 'anchor_ap_pref', 'anchor_dv_pref', ...
    'anchor_ap_weight', 'anchor_dv_weight', 'anchor_source', 'anchor_notes', ...
    'hierarchy_depth'});
end

function fig = draw_box_layout(labelsPlot, majorGroupPlot, centerXY, Aplot, opts, layoutInfo)
boxX = centerXY(:,1) - opts.boxWidth / 2;
boxY = centerXY(:,2) - opts.boxHeight / 2;

pixPerUnit = 70;
figW = max(1200, ceil(layoutInfo.totalWidth * pixPerUnit));
figH = max(900, ceil(layoutInfo.totalHeight * pixPerUnit));

fig = figure('Color', 'w', ...
    'Visible', opts.figureVisible, ...
    'Renderer', 'painters', ...
    'Units', 'pixels', ...
    'Position', [100 100 figW figH]);
ax = axes('Parent', fig);
set(ax, 'FontName', opts.fontName);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');

if opts.drawEdges
    [src, dst] = find(Aplot);
    if ~isempty(src)
        if numel(src) > opts.maxEdgesToDraw
            d2 = sum((centerXY(src,:) - centerXY(dst,:)).^2, 2);
            [~, order] = sort(d2, 'ascend');
            keepEdge = order(1:opts.maxEdgesToDraw);
            src = src(keepEdge);
            dst = dst(keepEdge);
        end

        for i = 1:numel(src)
            line(ax, [centerXY(src(i),1), centerXY(dst(i),1)], ...
                [centerXY(src(i),2), centerXY(dst(i),2)], ...
                'Color', opts.edgeColor, ...
                'LineWidth', opts.edgeLineWidth);
        end
    end
end

for i = 1:numel(labelsPlot)
    faceColor = resolve_box_face_color(majorGroupPlot{i}, opts);
    rectangle(ax, ...
        'Position', [boxX(i), boxY(i), opts.boxWidth, opts.boxHeight], ...
        'Curvature', opts.boxCurvature, ...
        'FaceColor', faceColor, ...
        'EdgeColor', opts.boxEdgeColor, ...
        'LineWidth', opts.lineWidth);

    localFontSize = opts.fontSize;
    if numel(labelsPlot{i}) >= 6
        localFontSize = max(4.5, opts.fontSize - 1);
    end

    text(ax, centerXY(i,1), centerXY(i,2), labelsPlot{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', localFontSize, ...
        'FontName', opts.fontName, ...
        'Interpreter', 'none');
end

padX = max(opts.boxWidth, 0.5);
padY = max(opts.boxHeight, 0.5);
xlim(ax, [min(boxX) - padX, max(boxX) + opts.boxWidth + padX]);
ylim(ax, [min(boxY) - padY, max(boxY) + opts.boxHeight + padY]);
end

function faceColor = resolve_box_face_color(majorGroup, opts)
if ~opts.colorByMajorGroup
    faceColor = opts.boxFaceColor;
    return;
end

switch majorGroup
    case 'cortical'
        faceColor = [0.90 0.94 1.00];
    case 'thalamic'
        faceColor = [1.00 0.92 0.80];
    case 'basal_ganglia'
        faceColor = [0.88 0.96 0.86];
    case 'diencephalon'
        faceColor = [0.98 0.89 0.94];
    case 'lower'
        faceColor = [0.92 0.92 0.92];
    otherwise
        faceColor = [1.00 1.00 1.00];
end
end

function files = export_outputs(fig, outPrefix, tbl, opts, cacheFile)
files = struct();
files.png = [outPrefix '.png'];
files.pdf = [outPrefix '.pdf'];
files.svg = '';
files.csv = [outPrefix '_coordinates.csv'];
files.mdsCache = cacheFile;

writetable(tbl, files.csv);
exportgraphics(fig, files.png, 'Resolution', opts.pngResolution);

try
    exportgraphics(fig, files.pdf, 'ContentType', 'vector');
catch
    print(fig, files.pdf, '-dpdf', '-painters');
end

if opts.exportSvg
    svgFile = [outPrefix '.svg'];
    try
        exportgraphics(fig, svgFile, 'ContentType', 'vector');
        files.svg = svgFile;
    catch
        try
            print(fig, svgFile, '-dsvg', '-painters');
            files.svg = svgFile;
        catch ME
            warning('make_modha_mds_box_layout:SvgExportFailed', ...
                'SVG export failed: %s', ME.message);
        end
    end
end
end
