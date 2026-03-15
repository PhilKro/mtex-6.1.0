function [ebsd] = pseudoSymmetryCorrection(ebsd, pseudoSym, varargin)
% PSEUDOSYMMETRYCORRECTION Corrects pseudo-symmetry artifacts using graph-based clustering.
%
% Strategy:
%   1. Calculate grains internally using a tight boundary setting.
%   2. Identify boundaries with misorientations matching the pseudo-symmetry.
%   3. Construct a graph and cluster connected grains.
%   4. Calculate metrics (size, ratio, pseudoSymBoundaryFraction) vectorized.
%   5. Identify "speckles" using a customizable function handle.
%   6. Select a "host" grain for each cluster using a customizable scoring function.
%   7. Rotate speckles to match the host orientation.
%   8. Clear temporary grain IDs and return the corrected EBSD data.
%
% Authors: Philipp Kroeker and Gemini (Google AI)
% Note: Collaboratively refactored for modularity, vectorization, and speed.
%
% Inputs:
%   ebsd      - @EBSD object
%   pseudoSym - @rotation (list of pseudo-symmetry operators)
%
% Optional Name-Value Pairs:
%   'SpeckleCondition'  - Function handle: @(metrics) returning logical array
%                         Example (purely size < 50px): @(m) m.size < 50 & m.pseudoSymBoundaryFraction > 0
%   'HostScore'         - Function handle: @(metrics) returning numeric score array
%   'CalcGrainAngle'    - Misorientation angle for initial grain calculation (default 5*degree)
%   'MisorientationTol' - Tolerance for pseudo-symmetry (default 5*degree)
%   'RotationThreshold' - Minimum misorientation to host to apply rotation (default 10*degree)
%   'UseMAD'            - Boolean to use MAD logic (lower MAD is identified as host) for hosts (default false)
%
% Outputs:
%   ebsd   - @EBSD object with corrected orientations and cleared grainId field

    %% 0. Parse Inputs & Setup Modularity
    p = inputParser;
    
    % Default condition for identifying speckles
    defaultSpeckleCond = @(m) (m.ratio > 0.1) & (m.pseudoSymBoundaryFraction > 0.3);
    
    % Default scoring for hosts (largest grain wins)
    defaultHostScore = @(m) m.size; 
    
    addRequired(p, 'ebsd');
    addRequired(p, 'pseudoSym');
    addParameter(p, 'SpeckleCondition', defaultSpeckleCond);
    addParameter(p, 'HostScore', defaultHostScore);
    addParameter(p, 'CalcGrainAngle', 5*degree);
    addParameter(p, 'MisorientationTol', 5*degree);
    addParameter(p, 'RotationThreshold', 10*degree);
    addParameter(p, 'UseMAD', false);
    
    parse(p, ebsd, pseudoSym, varargin{:});
    opts = p.Results;

    %% 1. Calculate grains with tight boundary for pseudo-symmetry identification
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', opts.CalcGrainAngle, 'boundary', 'tight');
    maxId = max(grains.id);

    %% 2. Identify Pseudo-Symmetry Boundaries & Build Graph
    gB = grains.boundary('indexed', 'indexed');
    
    isPseudo = false(size(gB));
    for k = 1:length(pseudoSym)
        isPseudo = isPseudo | (angle(gB.misorientation, pseudoSym(k)) < opts.MisorientationTol);
    end
    
    gB_ps = gB(isPseudo);
    
    if isempty(gB_ps)
        fprintf('No pseudo-symmetry boundaries found.\n');
        return;
    end
    
    % figure;
    % plot(gB)
    % hold on; plot(gB_ps, 'lineColor','r')

    edges = gB_ps.grainId;
    if ~isempty(edges)
        maxId = max(maxId, max(edges(:)));
    end
    
    % Build graph and extract cluster bins
    G = graph(edges(:,1), edges(:,2), [], maxId);
    bins = conncomp(G)'; % Transpose to column vector (maxId x 1)
    
    % rng(42);
    % cmap = hsv(50);
    % figure;
    % plot(grains, mod(bins - 1, 50) + 1, 'micronbar', 'off');
    % colormap(cmap(randperm(50), :));
    % drawnow;

    %% 3. Calculate Metrics & Identify Speckles
    metrics = struct();
    
    % A. Calculate Total Perimeter per Grain (Avoiding large array concatenation)
    all_gB = grains.boundary;
    ids_all = all_gB.grainId;
    len_all = all_gB.segLength;
    
    v1 = ids_all(:,1) > 0;
    v2 = ids_all(:,2) > 0;
    
    totalPerimeter = accumarray(ids_all(v1, 1), len_all(v1), [maxId, 1]) + ...
                     accumarray(ids_all(v2, 2), len_all(v2), [maxId, 1]);
   
    % B. Calculate Pseudo-Symmetry Perimeter per Grain
    ids_ps = gB_ps.grainId;
    len_ps = gB_ps.segLength;
    
    vps1 = ids_ps(:,1) > 0;
    vps2 = ids_ps(:,2) > 0;
    
    pseudoPerimeter = accumarray(ids_ps(vps1, 1), len_ps(vps1), [maxId, 1]) + ...
                      accumarray(ids_ps(vps2, 2), len_ps(vps2), [maxId, 1]);
    
    % C. Populate Metrics Struct
    metrics.size = zeros(maxId, 1);
    metrics.size(grains.id) = grains.numPixel;
    
    metrics.ratio = totalPerimeter ./ metrics.size;
    metrics.ratio(isinf(metrics.ratio) | isnan(metrics.ratio)) = 0;
    
    metrics.pseudoSymBoundaryFraction = pseudoPerimeter ./ totalPerimeter;
    metrics.pseudoSymBoundaryFraction(isnan(metrics.pseudoSymBoundaryFraction)) = 0; 
    
    % Pre-calculate the maximum size in each cluster for relative comparisons
    clusterMaxSize = accumarray(bins, metrics.size, [], @max);
    metrics.clusterMaxSize = clusterMaxSize(bins);
    
    % Handle MAD if requested
    if opts.UseMAD && (isfield(ebsd.prop, 'MAD') || isfield(ebsd.prop, 'mad'))
        madProp = 'MAD'; if isfield(ebsd.prop, 'mad'), madProp = 'mad'; end
        validID = ebsd.grainId > 0;
        
        % Calculate mean MAD per grain
        mad_max_id = max(maxId, max(ebsd.grainId(validID)));
        grainMeanMAD = accumarray(ebsd.grainId(validID), ebsd.prop.(madProp)(validID), [mad_max_id, 1], @mean, NaN);
        
        metrics.mad = grainMeanMAD(1:maxId);
        
        % If user didn't provide a custom HostScore, use the MAD logic:
        % Maximize score by penalizing MAD, constrained to top 20% of cluster size
        if isequal(opts.HostScore, defaultHostScore)
            opts.HostScore = @(m) -m.mad - 1e6 * (m.size < 0.2 * m.clusterMaxSize);
        end
    end

    % Evaluate the SpeckleCondition
    isSpeckle = opts.SpeckleCondition(metrics);
    
    % figure; plot(grains); hold on; plot(grains(isSpeckle), 'FaceColor','r')
    
    %% 4. Vectorized Host Selection
    grain_present_mask = false(maxId, 1);
    grain_present_mask(grains.id) = true;
    
    % Evaluate the modular HostScore
    scores = opts.HostScore(metrics);
    
    % Non-existent grains cannot be hosts.
    % (Removed the isSpeckle exclusion so a cluster composed entirely of speckles
    % still successfully selects a host based on the highest score).
    scores(~grain_present_mask) = -Inf;
    
    % rng(42);
    % cmap = hsv(50);
    % figure;
    % plot(grains, mod(scores - 1, 50) + 1, 'micronbar', 'off');
    % colormap(cmap(randperm(50), :));
    % drawnow;

    % Vectorized finding of the max score per bin using sortrows
    gIds = (1:maxId)';
    
    % Sort by bin (ascending) and then by score (descending)
    [~, sortIdx] = sortrows([bins, -scores]);
    sortedBins = bins(sortIdx);
    sortedIds  = gIds(sortIdx);
    
    % The first ID in each bin group has the highest score
    [uniqueBins, firstIdx] = unique(sortedBins, 'stable');
    
    % Filter out invalid bins (where all candidates had -Inf score)
    validBinsMask = (uniqueBins > 0) & (scores(sortedIds(firstIdx)) ~= -Inf);
    
    % Map the winning host ID back to the clusters
    hostForBin = zeros(max(bins), 1);
    actualHosts = sortedIds(firstIdx(validBinsMask));
    hostForBin(uniqueBins(validBinsMask)) = actualHosts;
    
    % Force the chosen hosts to NOT be speckles, protecting their boundaries
    % with other legitimate, non-speckle grains in the cluster.
    isSpeckle(actualHosts) = false;
    % figure; plot(grains); hold on; plot(grains(isSpeckle), 'FaceColor','r')
    
    % Map cluster hosts to individual grains
    host_assignments = hostForBin(bins);
    
    % Identify which speckles actually need mapping
    speckles_to_map_mask = isSpeckle & (host_assignments > 0);
    % figure; plot(grains); hold on; plot(grains(speckles_to_map_mask), 'FaceColor','b')

    %% 5. Vectorized Rotation Calculation
    rotations = orientation.nan(maxId, 1, pseudoSym(1).CS, pseudoSym(1).SS);
    
    grainOrientations = orientation.nan(maxId, 1, pseudoSym(1).CS);
    grainOrientations(grains.id) = grains.meanOrientation;
    
    s_ids = find(speckles_to_map_mask);
    
    if ~isempty(s_ids)
        h_ids = host_assignments(s_ids);
        
        ori_s = grainOrientations(s_ids);
        ori_h = grainOrientations(h_ids);
        
        % Only apply rotation if misorientation is significant
        needs_rot = angle(ori_s, ori_h) >= opts.RotationThreshold;
        
        s_ids_rot = s_ids(needs_rot);
        ori_s_rot = ori_s(needs_rot);
        ori_h_rot = ori_h(needs_rot);
        
        if ~isempty(s_ids_rot)
            mori = inv(ori_s_rot) .* ori_h_rot;
            
            all_syms = [pseudoSym, inv(pseudoSym)];
            dists = angle(mori, all_syms);
            [~, sym_idx] = min(dists, [], 2);
            
            rotations(s_ids_rot) = all_syms(sym_idx);
        end
    end
    
    %% 6. Update EBSD Data & Cleanup
    valid_grain_mask = ebsd.grainId > 0;
    pixel_grain_ids = ebsd.grainId(valid_grain_mask);
    
    pixel_rots = rotations(pixel_grain_ids);
    to_update = ~isnan(pixel_rots);
    
    if any(to_update)
        valid_indices = find(valid_grain_mask);
        update_indices = valid_indices(to_update);
        
        ebsd(update_indices).orientations = ebsd(update_indices).orientations .* pixel_rots(to_update);
        fprintf('Corrected %d pixels in pseudo-symmetry artifacts.\n', length(update_indices));
    end
    
    % Clear the ebsd.grainId field so the output object is clean
    if isfield(ebsd.prop, 'grainId')
        ebsd.prop = rmfield(ebsd.prop, 'grainId');
    end
end