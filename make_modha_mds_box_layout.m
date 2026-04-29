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
[centerXY, gridRC, layoutInfo] = assign_box_layout(Yplot, opts);

boxX = centerXY(:,1) - opts.boxWidth / 2;
boxY = centerXY(:,2) - opts.boxHeight / 2;

tbl = build_output_table(nodeIndex, labels, keepMask, Yplot, gridRC, centerXY, ...
    boxX, boxY, inDegree, outDegree, rootIndexAll, rootLabelAll, depthAll, majorGroupAll);

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
OUT.boxCenterXY = centerXY;
OUT.gridRC = gridRC;
OUT.parentOf = parentOfAll;
OUT.majorGroup = majorGroupAll;
OUT.inDegree = inDegree;
OUT.outDegree = outDegree;
OUT.stress = stress;
OUT.distanceInfo = distInfo;
OUT.mdsInfo = mdsInfo;
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

if ~strcmp(opts.layout, 'grid') && ~strcmp(opts.layout, 'circle')
    error('make_modha_mds_box_layout:BadLayout', ...
        'opts.layout must be ''grid'' or ''circle''.');
end

if ~strcmp(opts.mdsDisplay, 'off') && ~strcmp(opts.mdsDisplay, 'final') && ~strcmp(opts.mdsDisplay, 'iter')
    error('make_modha_mds_box_layout:BadMDSDisplay', ...
        'opts.mdsDisplay must be ''off'', ''final'', or ''iter''.');
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

function [centerXY, gridRC, layoutInfo] = assign_box_layout(Y, opts)
nNodes = size(Y, 1);
cellW = opts.boxWidth + opts.boxGapX;
cellH = opts.boxHeight + opts.boxGapY;

shapeList = candidate_shapes(nNodes, opts, cellW, cellH);
angles = linspace(0, pi, max(1, opts.nAngleSamples) + 1);
angles(end) = [];
if isempty(angles)
    angles = 0;
end

bestCost = inf;
bestCenterXY = [];
bestGridRC = [];
bestInfo = struct();

for s = 1:size(shapeList, 1)
    nRows = shapeList(s, 1);
    nCols = shapeList(s, 2);
    [candidateXY, candidateRC] = make_candidate_cells(nRows, nCols, nNodes, opts, cellW, cellH);

    for a = 1:numel(angles)
        Yrot = rotate_points(Y, angles(a));
        Yfit = fit_points_to_candidates(Yrot, candidateXY, opts.fitFill);
        [assignIdx, totalCost, method] = assign_points_to_cells(Yfit, candidateXY, opts);

        if totalCost < bestCost
            bestCost = totalCost;
            bestCenterXY = candidateXY(assignIdx, :);
            bestGridRC = candidateRC(assignIdx, :);
            bestInfo.nRows = nRows;
            bestInfo.nCols = nCols;
            bestInfo.angleDeg = angles(a) * 180 / pi;
            bestInfo.assignmentMethod = method;
            bestInfo.totalCost = totalCost;
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

function [assignIdx, totalCost, method] = assign_points_to_cells(Yfit, candidateXY, opts)
costMatrix = pdist2(Yfit, candidateXY, 'squaredeuclidean');

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
        totalCost = sum(costMatrix(sub2ind(size(costMatrix), (1:size(costMatrix,1))', assignIdx)));
        method = 'matchpairs';
        return;
    end
end

[assignIdx, totalCost] = greedy_assignment(costMatrix, Yfit);
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

function tbl = build_output_table(nodeIndex, labels, keepMask, Yplot, gridRC, centerXY, ...
    boxX, boxY, inDegree, outDegree, rootIndexAll, rootLabelAll, depthAll, majorGroupAll)
nNodes = numel(nodeIndex);

mdsX = nan(nNodes, 1);
mdsY = nan(nNodes, 1);
gridRow = nan(nNodes, 1);
gridCol = nan(nNodes, 1);
boxCenterX = nan(nNodes, 1);
boxCenterY = nan(nNodes, 1);
boxLeftX = nan(nNodes, 1);
boxBottomY = nan(nNodes, 1);

plotIdx = find(keepMask);
mdsX(plotIdx) = Yplot(:,1);
mdsY(plotIdx) = Yplot(:,2);
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
    depthAll(:), ...
    'VariableNames', { ...
    'index', 'acronym', 'is_plotted', 'mds_x', 'mds_y', 'grid_row', 'grid_col', ...
    'box_center_x', 'box_center_y', 'box_x', 'box_y', ...
    'in_degree', 'out_degree', 'root_index', 'root_acronym', 'major_group', 'hierarchy_depth'});
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
