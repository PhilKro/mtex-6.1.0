function [grains, ebsd] = graphedPseudoSymRemoval(ebsd, grains, pseudoSym, ratioThreshold, fractionThreshold, varargin)
% GRAPHEDPSEUDOSYMREMOVAL Removes EBSD artifacts using a graph-based perimeter fraction approach.
%
% New Strategy:
%   1. Build a graph connecting all grains with pseudo-symmetry boundaries.
%   2. Find connected components (clusters) in this graph.
%   3. Within each cluster, identify "speckle" grains based on thresholds.
%   4. Merge speckles into a designated "host" grain within the same cluster.
%   5. The host is the largest non-speckle grain, or the "best" speckle 
%      (lowest MAD/largest size) if all are speckles.
%
% Inputs:
%   ebsd              - EBSD data
%   grains            - grain2d object
%   pseudoSym         - list of pseudo-symmetry rotations
%   ratioThreshold    - (Optional) Threshold for Perimeter/Size ratio (default 0.37)
%   fractionThreshold - (Optional) Threshold for PseudoBoundary/TotalBoundary (default 0.5)
%   'disregardMAD'    - (Flag) Disregard MAD values for host determination
%
% Outputs:
%   grains - Updated grains
%   ebsd   - Updated EBSD data (orientations rotated)

    if nargin < 4 || isempty(ratioThreshold)
        ratioThreshold = 0.2;
    end
    if nargin < 5 || isempty(fractionThreshold)
        fractionThreshold = 0.3;
    end
    
    disregardMAD = any(strcmpi(varargin, 'disregardMAD'));

    maxId = max(grains.id);

    %% 1. Identify Pseudosymmetry Boundaries & Build Graph
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
    bins = conncomp(G); % `bins(i)` is the cluster ID for grain with ID `i`

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
    %% 2. Calculate Metrics & Identify Speckles
    % A. Total Perimeter per Grain
    all_gB = grains.boundary;
    ids_all = all_gB.grainId;
    len_all = all_gB.segLength;
    
    valid_mask = ids_all > 0;
    flat_ids = [ids_all(valid_mask(:,1), 1); ids_all(valid_mask(:,2), 2)];
    flat_len = [len_all(valid_mask(:,1));    len_all(valid_mask(:,2))];
    
    totalPerimeter = accumarray(flat_ids, flat_len, [maxId, 1]);
    
    % B. Pseudosymmetry Perimeter per Grain
    ids_ps = gB_ps.grainId;
    len_ps = gB_ps.segLength;
    
    flat_ids_ps = [ids_ps(ids_ps(:,1)>0,1); ids_ps(ids_ps(:,2)>0,2)];
    flat_len_ps = [len_ps(ids_ps(:,1)>0);    len_ps(ids_ps(:,2)>0)];

    pseudoPerimeter = accumarray(flat_ids_ps, flat_len_ps, [maxId, 1]);
    
    % C. Calculate Ratios and Identify Speckles
    g_numPixel = zeros(maxId, 1);
    g_numPixel(grains.id) = grains.numPixel;

    g_ratio = totalPerimeter ./ g_numPixel;
    g_ratio(isinf(g_ratio) | isnan(g_ratio)) = 0;
    
    g_fraction = pseudoPerimeter ./ totalPerimeter;
    g_fraction(isnan(g_fraction)) = 0; % Handle grains with 0 perimeter (singularities)
    
    isSpeckle = (g_ratio > ratioThreshold) & (g_fraction > fractionThreshold);

    %% 3. Process Clusters: Determine Rotations and Merges
    rotations = orientation.nan(maxId, 1, pseudoSym(1).CS, pseudoSym(1).SS);
    gB_to_merge_mask = false(size(gB_ps));

    % Check for MAD property to use as a tie-breaker
    useMAD = false;
    if ~disregardMAD && (isfield(ebsd.prop, 'MAD') || isfield(ebsd.prop, 'mad'))
        warning('MAD implementation is shakey!')
        useMAD = true;
        madProp = 'MAD'; if isfield(ebsd.prop, 'mad'), madProp = 'mad'; end
        validID = ebsd.grainId > 0;
        mad_max_id = max(maxId, max(ebsd.grainId(validID)));
        grainMeanMAD = accumarray(ebsd.grainId(validID), ebsd.prop.(madProp)(validID), [mad_max_id, 1], @mean, NaN);
    end

    unique_bins = unique(bins(bins>0));
    
    for i = 1:length(unique_bins)
        b = unique_bins(i);
        cluster_ids = find(bins == b);
        cluster_ids = cluster_ids(ismember(cluster_ids, grains.id));
        
        if length(cluster_ids) < 2, continue; end
        
        % Identify speckles and hosts within the cluster
        is_speckle_in_cluster = isSpeckle(cluster_ids);
        speckle_ids = cluster_ids(is_speckle_in_cluster);
        
        if isempty(speckle_ids), continue; end % No speckles to remove in this cluster

        % Determine the single host for the cluster
        host_id = 0;
        candidates = cluster_ids;
        
        if ~isempty(candidates)
            cand_sizes = g_numPixel(candidates);
            
            if useMAD
                % MAD is main criterion (if not disregarded)
                % Constraint: Host must not be "comparatively small" (Top ~20% of candidates)
                max_s = max(cand_sizes);
                is_large_enough = cand_sizes >= 0.2 * max_s;
                
                valid_cands = candidates(is_large_enough);
                [~, min_mad_idx] = min(grainMeanMAD(valid_cands));
                host_id = valid_cands(min_mad_idx);
            else
                % Size is main criterion
                [~, max_idx] = max(cand_sizes);
                host_id = candidates(max_idx);
            end
        end
        
        if host_id == 0, continue; end

        % Get host orientation
        ori_host = grains(grains.id2ind(host_id)).meanOrientation;
        
        % For troubleshooting
        % fig = figure;
        % plot(grains(grains.id2ind(cluster_ids)), grains(grains.id2ind(cluster_ids)).meanOrientation)
        % hold on
        % nextAxis
        % plot(grains(grains.id2ind(host_id)), 'FaceColor','RebeccaPurple')
        % hold on
        % plot(grains(grains.id2ind(speckle_ids)), 'FaceColor','b')
        % nextAxis
        % plot(grains(grains.id2ind(cluster_ids)))
        % hold on
        % plot(grains(grains.id2ind(host_id)), 'FaceColor','RebeccaPurple')

        % Determine rotations for all other grains to be merged into the host
        all_syms = [pseudoSym, inv(pseudoSym)]; % Consider inverse symmetries too
        for k = 1:length(cluster_ids)
            s_id = cluster_ids(k);
            if s_id == host_id, continue; end
            
            % Only rotate and merge if the grain is actually a speckle
            if ~isSpeckle(s_id)
                continue;
            end
            
            ori_speckle = grains(grains.id2ind(s_id)).meanOrientation;

           % Exclude grains that already have an orientation close to the host
            if angle(ori_speckle, ori_host) < 10*degree
                hold on
                % For troubleshooting
                % plot(grains(grains.id2ind(s_id)), 'FaceColor','w')
                continue;
            end
            
            % Find optimal symmetry to match host
            mori = inv(ori_speckle) * ori_host;
            
            [~, sym_idx] = min(angle(mori, all_syms));
            rotations(s_id) = all_syms(sym_idx);
            hold on

            % For troubleshooting
            % plot(grains(grains.id2ind(s_id)), 'FaceColor','k')
        end
        % For troubleshooting
        % close(fig)

        % Flag all internal pseudo-sym boundaries within this cluster for merging
        % Only merge boundaries where at least one side is a speckle.
        % (Preserves boundaries between two "real" grains even if they are pseudo-symmetric)
        ids_in_cluster = ismember(gB_ps.grainId(:,1), cluster_ids) & ismember(gB_ps.grainId(:,2), cluster_ids);
        at_least_one_speckle = isSpeckle(gB_ps.grainId(:,1)) | isSpeckle(gB_ps.grainId(:,2));
        
        gB_cluster_mask = ids_in_cluster & at_least_one_speckle;
        gB_to_merge_mask = gB_to_merge_mask | gB_cluster_mask;
    end

    %% 4. Update EBSD data and Merge Grains
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
end