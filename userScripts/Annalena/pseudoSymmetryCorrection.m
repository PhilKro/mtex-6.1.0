function [ebsd] = pseudoSymmetryCorrection(ebsd, pseudoSym, varargin)
% PSEUDOSYMMETRYCORRECTION Corrects pseudo-symmetry artifacts using graph-based clustering.
%
% Available Metrics for Function Handles (m.XXX):
%   m.size                      : Number of pixels in the grain.
%   m.ratio                     : Total boundary length / size.
%   m.pseudoSymBoundaryFraction : Pseudo-symmetry boundary length / Total boundary length.
%   m.tortuosity                : Pseudo-symmetry boundary length / Bounding box diagonal.
%   m.clusterMaxSize            : Size of the largest grain in the same pseudo-symmetry cluster.
%   m.mad                       : Mean Angular Deviation (requires ebsd.prop.mad).
%
% Strategy:
%   1. Calculate grains internally using a tight boundary setting.
%   2. Identify boundaries with misorientations matching the pseudo-symmetry.
%   3. Construct a graph and cluster connected grains.
%   4. Dynamically evaluate and calculate only the requested metrics.
%   5. Identify "speckles" using a customizable function handle.
%   6. Select a "host" grain for each cluster using a customizable scoring function.
%   7. Rotate speckles to match the host orientation.
%   8. Clear temporary grain IDs and return the corrected EBSD data.
%
% Authors: Philipp Kroeker and Gemini (Google AI)
%
% Inputs:
%   ebsd      - @EBSD object
%   pseudoSym - @rotation (list of pseudo-symmetry operators)
%
% Optional Name-Value Pairs:
%   'SpeckleCondition'  - Function handle: @(metrics) returning logical array
%                         Example (purely size < 50px): @(m) m.size < 50 & m.pseudoSymBoundaryFraction > 0
%                         Example (tortuosity): @(m) m.tortuosity > 1.5 & m.pseudoSymBoundaryFraction > 0.5
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
        if isfield(ebsd.prop, 'grainId')
            ebsd.prop = rmfield(ebsd.prop, 'grainId');
        end
        return;
    end

    edges = gB_ps.grainId;
    if ~isempty(edges)
        maxId = max(maxId, max(edges(:)));
    end
    
    % Build graph and extract cluster bins
    G = graph(edges(:,1), edges(:,2), [], maxId);
    bins = conncomp(G)'; % Transpose to column vector (maxId x 1)

    %% 3. Dynamically Calculate Metrics
    % If UseMAD is requested and no custom HostScore is provided, inject the MAD logic
    if opts.UseMAD && isequal(opts.HostScore, defaultHostScore)
        opts.HostScore = @(m) -m.mad - 1e6 * (m.size < 0.2 * m.clusterMaxSize);
    end

    % Inspect function strings to determine required calculations
    condStr = [func2str(opts.SpeckleCondition), ' ', func2str(opts.HostScore)];
    
    needSize        = contains(condStr, '.size');
    needRatio       = contains(condStr, '.ratio');
    needFraction    = contains(condStr, '.pseudoSymBoundaryFraction');
    needTortuosity  = contains(condStr, '.tortuosity');
    needClusterMax  = contains(condStr, '.clusterMaxSize');
    needMad         = contains(condStr, '.mad');
    
    % Resolve calculation dependencies
    needTotalPerim  = needRatio || needFraction;
    needPseudoPerim = needFraction || needTortuosity;
    needSize        = needSize || needRatio || needClusterMax;
    
    metrics = struct();
    
    % A. Grain Size
    if needSize
        metrics.size = zeros(maxId, 1);
        metrics.size(grains.id) = grains.numPixel;
    end
    
    % B. Total Perimeter & Ratio
    if needTotalPerim
        all_gB = grains.boundary;
        ids_all = all_gB.grainId;
        len_all = all_gB.segLength;
        
        v1 = ids_all(:,1) > 0;
        v2 = ids_all(:,2) > 0;
        
        totalPerimeter = accumarray(ids_all(v1, 1), len_all(v1), [maxId, 1]) + ...
                         accumarray(ids_all(v2, 2), len_all(v2), [maxId, 1]);
                         
        if needRatio
            metrics.ratio = totalPerimeter ./ metrics.size;
            metrics.ratio(isinf(metrics.ratio) | isnan(metrics.ratio)) = 0;
        end
    end
    
    % C. Pseudo-Symmetry Perimeter, Fraction & Tortuosity
    if needPseudoPerim
        ids_ps = gB_ps.grainId;
        len_ps = gB_ps.segLength;
        
        vps1 = ids_ps(:,1) > 0;
        vps2 = ids_ps(:,2) > 0;
        
        flat_ids_ps = [ids_ps(vps1, 1); ids_ps(vps2, 2)];
        flat_len_ps = [len_ps(vps1); len_ps(vps2)];
        
        pseudoPerimeter = accumarray(flat_ids_ps, flat_len_ps, [maxId, 1]);
        
        if needFraction
            metrics.pseudoSymBoundaryFraction = pseudoPerimeter ./ totalPerimeter;
            metrics.pseudoSymBoundaryFraction(isnan(metrics.pseudoSymBoundaryFraction)) = 0; 
        end
        
        if needTortuosity
            mid_x_ps = gB_ps.midPoint.x;
            mid_y_ps = gB_ps.midPoint.y;
            flat_x_ps = [mid_x_ps(vps1); mid_x_ps(vps2)];
            flat_y_ps = [mid_y_ps(vps1); mid_y_ps(vps2)];

            min_x = accumarray(flat_ids_ps, flat_x_ps, [maxId, 1], @min, NaN);
            max_x = accumarray(flat_ids_ps, flat_x_ps, [maxId, 1], @max, NaN);
            min_y = accumarray(flat_ids_ps, flat_y_ps, [maxId, 1], @min, NaN);
            max_y = accumarray(flat_ids_ps, flat_y_ps, [maxId, 1], @max, NaN);
            
            diag_length = sqrt((max_x - min_x).^2 + (max_y - min_y).^2);
            
            metrics.tortuosity = pseudoPerimeter ./ diag_length;
            metrics.tortuosity(isnan(metrics.tortuosity) | isinf(metrics.tortuosity)) = 0;
        end
    end
    
    % D. Cluster Maximum Size
    if needClusterMax
        clusterMaxSize = accumarray(bins, metrics.size, [], @max);
        metrics.clusterMaxSize = clusterMaxSize(bins);
    end
    
    % E. Mean Angular Deviation (MAD)
    if needMad
        madProp = 'MAD'; if isfield(ebsd.prop, 'mad'), madProp = 'mad'; end
        if isfield(ebsd.prop, madProp)
            validID = ebsd.grainId > 0;
            mad_max_id = max(maxId, max(ebsd.grainId(validID)));
            grainMeanMAD = accumarray(ebsd.grainId(validID), ebsd.prop.(madProp)(validID), [mad_max_id, 1], @mean, NaN);
            metrics.mad = grainMeanMAD(1:maxId);
        else
            warning('MAD property not found in EBSD data. Defaulting metrics.mad to NaN.');
            metrics.mad = nan(maxId, 1);
        end
    end

    % Evaluate the SpeckleCondition
    isSpeckle = opts.SpeckleCondition(metrics);
    
    %% 4. Vectorized Host Selection
    grain_present_mask = false(maxId, 1);
    grain_present_mask(grains.id) = true;
    
    scores = opts.HostScore(metrics);
    scores(~grain_present_mask) = -Inf;

    gIds = (1:maxId)';
    
    [~, sortIdx] = sortrows([bins, -scores]);
    sortedBins = bins(sortIdx);
    sortedIds  = gIds(sortIdx);
    
    [uniqueBins, firstIdx] = unique(sortedBins, 'stable');
    
    validBinsMask = (uniqueBins > 0) & (scores(sortedIds(firstIdx)) ~= -Inf);
    
    hostForBin = zeros(max(bins), 1);
    actualHosts = sortedIds(firstIdx(validBinsMask));
    hostForBin(uniqueBins(validBinsMask)) = actualHosts;
    
    isSpeckle(actualHosts) = false;
    
    host_assignments = hostForBin(bins);
    speckles_to_map_mask = isSpeckle & (host_assignments > 0);

    %% 5. Vectorized Rotation Calculation
    rotations = orientation.nan(maxId, 1, pseudoSym(1).CS, pseudoSym(1).SS);
    
    grainOrientations = orientation.nan(maxId, 1, pseudoSym(1).CS);
    grainOrientations(grains.id) = grains.meanOrientation;
    
    s_ids = find(speckles_to_map_mask);
    
    if ~isempty(s_ids)
        h_ids = host_assignments(s_ids);
        
        ori_s = grainOrientations(s_ids);
        ori_h = grainOrientations(h_ids);
        
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
    
    if isfield(ebsd.prop, 'grainId')
        ebsd.prop = rmfield(ebsd.prop, 'grainId');
    end
end