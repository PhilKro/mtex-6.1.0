function [grains, ebsd] = graphedPseudoSymRemoval(ebsd, grains, pseudoSym, ratioThreshold, fractionThreshold, varargin)
% GRAPHEDPSEUDOSYMREMOVAL Corrects pseudo-symmetry artifacts using graph-based clustering.
%
% Strategy:
%   1. Identify boundaries with misorientations matching the pseudo-symmetry.
%   2. Construct a graph where nodes are grains and edges are pseudo-symmetry boundaries.
%   3. Cluster connected grains (components).
%   4. Identify "speckles" (artifacts) within clusters using perimeter/area and 
%      pseudo-symmetry boundary fraction thresholds.
%   5. Select a "host" grain for each cluster (lowest MAD or largest size).
%   6. Rotate speckles to match the host orientation.
%   7. Merge speckles into the host.
%
% Inputs:
%   ebsd              - @EBSD object
%   grains            - @grain2d object
%   pseudoSym         - @rotation (list of pseudo-symmetry operators)
%   ratioThreshold    - (Optional) Threshold for Perimeter/Size ratio (default 0.37)
%   fractionThreshold - (Optional) Threshold for PseudoBoundary/TotalBoundary (default 0.5)
%   varargin          - 'disregardMAD': flag to ignore MAD for host selection.
%
% Outputs:
%   grains - @grain2d object (Note: geometry may be inconsistent, see warning)
%   ebsd   - @EBSD object with corrected orientations

    if nargin < 4 || isempty(ratioThreshold)
        ratioThreshold = 0.1;
    end
    if nargin < 5 || isempty(fractionThreshold)
        fractionThreshold = 0.3;
    end
    
    disregardMAD = any(strcmpi(varargin, 'disregardMAD'));

    maxId = max(grains.id);

    %% 1. Identify Pseudo-Symmetry Boundaries & Build Graph
    gB = grains.boundary('indexed', 'indexed');
    
    isPseudo = false(size(gB));
    for k = 1:length(pseudoSym)
        isPseudo = isPseudo | (angle(gB.misorientation, pseudoSym(k)) < 5*degree);
    end
    
    gB_ps = gB(isPseudo);
    
    if isempty(gB_ps)
        fprintf('No pseudo-symmetry boundaries found.\n');
        return;
    end
    
    edges = gB_ps.grainId;
    % Ensure maxId is large enough for graph nodes
    if ~isempty(edges)
        maxId = max(maxId, max(edges(:)));
    end

    G = graph(edges(:,1), edges(:,2), [], maxId);
    bins = conncomp(G);

    %% For testing: Plot the clusters
    % clusterIds = zeros(length(grains),1);
    % valid_ids_for_plot = grains.id(grains.id <= length(bins));
    % if ~isempty(valid_ids_for_plot)
    %     clusterIds(grains.id2ind(valid_ids_for_plot)) = bins(valid_ids_for_plot);
    % end
    % 
    % rng(42);
    % cmap = hsv(50);
    % 
    % figure;
    % plot(grains, mod(clusterIds - 1, 50) + 1, 'micronbar', 'off');
    % colormap(cmap(randperm(50), :));
    % title('Pseudo-Symmetry Clusters');
    % drawnow;
    %% 2. Calculate Metrics & Identify Speckles (Artifacts)
    % A. Calculate Total Perimeter per Grain
    all_gB = grains.boundary;
    ids_all = all_gB.grainId;
    len_all = all_gB.segLength;
    
    valid_mask = ids_all > 0;
    flat_ids = [ids_all(valid_mask(:,1), 1); ids_all(valid_mask(:,2), 2)];
    flat_len = [len_all(valid_mask(:,1));    len_all(valid_mask(:,2))];
    
    totalPerimeter = accumarray(flat_ids, flat_len, [maxId, 1]);
    
    % B. Calculate Pseudo-Symmetry Perimeter per Grain
    ids_ps = gB_ps.grainId;
    len_ps = gB_ps.segLength;
    
    flat_ids_ps = [ids_ps(ids_ps(:,1)>0,1); ids_ps(ids_ps(:,2)>0,2)];
    flat_len_ps = [len_ps(ids_ps(:,1)>0);    len_ps(ids_ps(:,2)>0)];

    pseudoPerimeter = accumarray(flat_ids_ps, flat_len_ps, [maxId, 1]);
    
    % C. Apply Thresholds to Identify Speckles
    g_numPixel = zeros(maxId, 1);
    g_numPixel(grains.id) = grains.numPixel;

    g_ratio = totalPerimeter ./ g_numPixel;
    g_ratio(isinf(g_ratio) | isnan(g_ratio)) = 0;
    
    g_fraction = pseudoPerimeter ./ totalPerimeter;
    g_fraction(isnan(g_fraction)) = 0; % Handle grains with 0 perimeter (singularities)
    
    isSpeckle = (g_ratio > ratioThreshold) & (g_fraction > fractionThreshold);

    %% 3. Process Clusters: Select Hosts & Assign Merges
    rotations = orientation.nan(maxId, 1, pseudoSym(1).CS, pseudoSym(1).SS);
    
    % Optimization: Pre-fetch orientations for O(1) lookup
    grainOrientations = orientation.nan(maxId, 1, pseudoSym(1).CS);
    grainOrientations(grains.id) = grains.meanOrientation;
    
    % Prepare MAD data for host selection (tie-breaker)
    useMAD = false;
    if ~disregardMAD && (isfield(ebsd.prop, 'MAD') || isfield(ebsd.prop, 'mad'))
        warning('MAD implementation is shakey!')
        useMAD = true;
        madProp = 'MAD'; if isfield(ebsd.prop, 'mad'), madProp = 'mad'; end
        validID = ebsd.grainId > 0;
        mad_max_id = max(maxId, max(ebsd.grainId(validID)));
        grainMeanMAD = accumarray(ebsd.grainId(validID), ebsd.prop.(madProp)(validID), [mad_max_id, 1], @mean, NaN);
    end

    % Optimization: Sort bins to group clusters efficiently (avoiding find())
    involved_nodes = unique(edges(:));
    relevant_bins_per_node = bins(involved_nodes);
    
    [sorted_bins, sort_idx] = sort(relevant_bins_per_node);
    sorted_nodes = involved_nodes(sort_idx);
    
    [~, i_start] = unique(sorted_bins, 'first');
    [~, i_end] = unique(sorted_bins, 'last');
    
    grain_present_mask = false(maxId, 1);
    grain_present_mask(grains.id) = true;
    
    host_assignments = zeros(maxId, 1);

    for i = 1:length(i_start)
        cluster_ids = sorted_nodes(i_start(i):i_end(i));
        cluster_ids = cluster_ids(grain_present_mask(cluster_ids));
        
        if length(cluster_ids) < 2, continue; end
        
        % Identify speckle candidates within the cluster
        is_speckle_in_cluster = isSpeckle(cluster_ids);
        speckle_ids = cluster_ids(is_speckle_in_cluster);
        
        if isempty(speckle_ids), continue; end

        % Determine the single host for the cluster
        host_id = 0;
        candidates = cluster_ids;
        
        if ~isempty(candidates)
            cand_sizes = g_numPixel(candidates);
            
            if useMAD
                % Criteria: Lowest MAD, constrained to top 20% size
                max_s = max(cand_sizes);
                is_large_enough = cand_sizes >= 0.2 * max_s;
                
                valid_cands = candidates(is_large_enough);
                [~, min_mad_idx] = min(grainMeanMAD(valid_cands));
                host_id = valid_cands(min_mad_idx);
            else
                % Criteria: Largest Size
                [~, max_idx] = max(cand_sizes);
                host_id = candidates(max_idx);
            end
        end
        
        if host_id == 0, continue; end

        % Mark speckles for merging into the selected host
        speckles_to_map = cluster_ids(cluster_ids ~= host_id & isSpeckle(cluster_ids));
        host_assignments(speckles_to_map) = host_id;
    end
    
    %% 4. Vectorized Rotation Calculation
    s_ids = find(host_assignments > 0);
    
    if ~isempty(s_ids)
        h_ids = host_assignments(s_ids);
        
        ori_s = grainOrientations(s_ids);
        ori_h = grainOrientations(h_ids);
        
        % Only apply rotation if misorientation is significant (>10 deg)
        needs_rot = angle(ori_s, ori_h) >= 10*degree;
        
        s_ids_rot = s_ids(needs_rot);
        ori_s_rot = ori_s(needs_rot);
        ori_h_rot = ori_h(needs_rot);
        
        if ~isempty(s_ids_rot)
            % Find best symmetry operator to match host orientation
            mori = inv(ori_s_rot) .* ori_h_rot;
            
            all_syms = [pseudoSym, inv(pseudoSym)];
            dists = angle(mori, all_syms);
            [~, sym_idx] = min(dists, [], 2);
            
            rotations(s_ids_rot) = all_syms(sym_idx);
        end
    end
    
    % Identify boundaries where at least one side is a speckle
    gB_to_merge_mask = isSpeckle(gB_ps.grainId(:,1)) | isSpeckle(gB_ps.grainId(:,2));

    %% 5. Update EBSD Data and Merge
    % A. Apply rotations to EBSD data
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
    
    % B. Merge Grains
    gB_to_merge = gB_ps(gB_to_merge_mask);
    
    if ~isempty(gB_to_merge)
        [grains, parentIdMap] = merge(grains, gB_to_merge);
        
        % Update ebsd.grainId using the map from old grain IDs to new grain IDs
        % This is safe as long as max(ebsd.grainId) is not larger than the map length
        ebsd.grainId(ebsd.grainId > 0) = parentIdMap(ebsd.grainId(ebsd.grainId > 0));
        fprintf('Merged %d pseudo-symmetry boundaries.\n', length(gB_to_merge));
    end

    % Final Warning
    warning('graphedPseudoSymRemoval:garbageGrains', ...
        'The grains output of this function is not fully consistent (garbage). Please recalculate grains using the returned ebsd variable.');
end