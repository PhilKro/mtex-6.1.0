% c:\Users\phkr\OneDrive - empa.ch\mtex-6.1.0\userScripts\Annalena\removePseudoSymmetry.m
function [grains, ebsd] = removePseudoSymmetry(ebsd, grains, pseudoSym, varargin)
% removePseudoSymmetry merges grains that are likely pseudo-symmetry artifacts
%
% Criteria for removal (AND/Weighted):
% 1. Boundary misorientation matches pseudoSym.
% 2. Boundary curvature is high (indicating inclusion/artifact).
% 3. Grain size is small.
%
% Options:
%  'curvatureWeight' - weight for curvature (default 1)
%  'sizeWeight'      - weight for inverse grain size (default 5)
%  'threshold'       - score threshold for merging (default 2.5)
%  'sizeOnly'        - bool, if true, decision is based only on size (default false)
%  'curvatureOnly'   - bool, if true, decision is based only on curvature (default false)
%  'ignoreMAD'       - bool, if true, disregard MAD values even if present (default false)

    % Default weights
    w_curv = 1; 
    w_size = 5; 
    threshold = 2.5;
    sizeOnly = false;
    curvatureOnly = false;
    ignoreMAD = false;

    % Parse options
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'curvatureWeight', w_curv = varargin{i+1};
            case 'sizeWeight', w_size = varargin{i+1};
            case 'threshold', threshold = varargin{i+1};
            case 'sizeOnly', sizeOnly = varargin{i+1};
            case 'curvatureOnly', curvatureOnly = varargin{i+1};
            case 'ignoreMAD', ignoreMAD = varargin{i+1};
        end
    end

    % 1. Identify Pseudo-Symmetry Boundaries
    % Use a tolerance (e.g., 5 degrees)
    gB = grains.boundary('indexed','indexed');
    isPseudo = false(size(gB));
    for k = 1:length(pseudoSym)
        isPseudo = isPseudo | (angle(gB.misorientation, pseudoSym(k)) < 5*degree);
    end
    gB_ps = gB(isPseudo);

    if isempty(gB_ps)
        return;
    end

    % Check for MAD property to determine which grain to overwrite
    % Calculate mean MAD per grain only if boundaries are found (Optimization)
    useMAD = false;
    madProp = '';
    if ~ignoreMAD
        if isfield(ebsd.prop, 'MAD')
            useMAD = true; madProp = 'MAD';
        elseif isfield(ebsd.prop, 'mad')
            useMAD = true; madProp = 'mad';
        end
    end
    
    grainMeanMAD = [];
    if useMAD
        validID = ebsd.grainId > 0;
        grainMeanMAD = accumarray(ebsd.grainId(validID), ebsd.prop.(madProp)(validID), [length(grains), 1], @mean, NaN);
    end
    
    % Get IDs of grains on both sides of the boundary
    grainIds = gB_ps.grainId;
    
    % Filter out outer boundaries (ID=0)
    validMask = all(grainIds > 0, 2);
    gB_ps = gB_ps(validMask);
    grainIds = grainIds(validMask, :);
    
    % Sort grain IDs to ensure consistent pairs (Grain A - Grain B)
    grainIds = sort(grainIds, 2);
    
    % Find unique grain pairs to aggregate curvature properties
    [uPairs, ~, pairIdx] = unique(grainIds, 'rows');

    % 2. Calculate Curvature (Aggregated per pair)
    try
        % Calculate curvature for all segments
        kappa_seg = abs(gB_ps.curvature);
        
        % Disregard curvature values < 1e-5 (essentially 0/flat)
        isValidCurv = kappa_seg > 1e-5;
        
        % Compute mean curvature per grain pair using only valid segments
        sumKappa = accumarray(pairIdx(isValidCurv), kappa_seg(isValidCurv), [size(uPairs,1), 1]);
        countKappa = accumarray(pairIdx(isValidCurv), 1, [size(uPairs,1), 1]);
        
        meanKappa = zeros(size(uPairs,1), 1);
        hasData = countKappa > 0;
        meanKappa(hasData) = sumKappa(hasData) ./ countKappa(hasData);
        
    catch
        warning('Curvature calculation failed. Ensure grains are smoothed.');
        meanKappa = zeros(size(uPairs,1), 1);
    end

    % 3. Grain Sizes (for the unique pairs)
    % Retrieve sizes (numPixel)
    % Map IDs to indices safely
    ind1 = grains.id2ind(uPairs(:,1));
    ind2 = grains.id2ind(uPairs(:,2));
    
    s1 = reshape(grains.numPixel(ind1), [], 1);
    s2 = reshape(grains.numPixel(ind2), [], 1);
    
    % We are interested in the size of the *smaller* grain (the potential inclusion)
    [minSize, minIdx] = min([s1, s2], [], 2);

    % 4. Calculate Score
    % Score increases with high curvature and small size (1/size)
    % Logic:
    % - High Curvature + Small Grain = High Score -> Merge (Pseudo-Symmetry Artifact)
    % - Low Curvature + Big Grain    = Low Score  -> Keep  (Real Twin)
    if sizeOnly
        
        w_curv = 0;
        fprintf('Note: Pseudo-symmetry cleanup running in size-only mode.\n');
    end
    if curvatureOnly
        w_size = 0;
        fprintf('Note: Pseudo-symmetry cleanup running in curvature-only mode.\n');
    end

    score = w_curv * meanKappa + w_size ./ minSize;
    
    % Identify boundaries to merge
    toMergePairs = score > threshold;
    
    % Map back to the original boundary segments
    pairsToMerge = uPairs(toMergePairs, :);
    
    % Select all segments that belong to the pairs identified for merging
    maskMerge = ismember(grainIds, pairsToMerge, 'rows');
    gB_to_merge = gB_ps(maskMerge);

    if isempty(gB_to_merge)
        return;
    end

    % 5. Update EBSD Orientations
    % Before merging, we want to rotate the 'bad' grains to match the parent.
    % We identify the bad grain as the one with higher MAD (if available) or smaller size.
    
    % Determine bad/good based on MAD or Size for each merging pair
    badIds = [];
    goodIds = [];
    
    for k = 1:size(pairsToMerge, 1)
        id1 = pairsToMerge(k, 1);
        id2 = pairsToMerge(k, 2);
        
        is1Bad = false;
        
        if useMAD && ~isnan(grainMeanMAD(id1)) && ~isnan(grainMeanMAD(id2))
            if grainMeanMAD(id1) > grainMeanMAD(id2)
                is1Bad = true;
            elseif grainMeanMAD(id1) < grainMeanMAD(id2)
                is1Bad = false;
            else
                % Tie-break with size (smaller is bad)
                is1Bad = grains(id1).numPixel < grains(id2).numPixel;
            end
        else
            % Size only (smaller is bad)
            is1Bad = grains(id1).numPixel < grains(id2).numPixel;
        end
        
        if is1Bad
            badIds = [badIds; id1];
            goodIds = [goodIds; id2];
        else
            badIds = [badIds; id2];
            goodIds = [goodIds; id1];
        end
    end
    
    % Unique bad grains to process
    [uBadIds, idx] = unique(badIds);
    uGoodIds = goodIds(idx);
    
    % Vectorize rotation to avoid slow loop over grains
    countChanged = 0;
    % Initialize with NaNs to indicate no rotation needed
    grainRotations = orientation.nan(length(grains), 1, pseudoSym(1).CS);
    
    for k = 1:length(uBadIds)
        bid = uBadIds(k);
        gid = uGoodIds(k);
        
        oriBad = grains(bid).meanOrientation;
        oriGood = grains(gid).meanOrientation;
        
        % Calculate required rotation: oriBad * rot ~ oriGood
        mori = inv(oriBad) * oriGood;
        
        % Find best matching symmetry from list
        bestRot = orientation.nan;
        minAng = inf;
        
        for s = 1:length(pseudoSym)
            sym = pseudoSym(s);
            
            d = angle(mori, sym);
            if d < minAng
                minAng = d;
                bestRot = sym;
            end
            
            % Check inverse (important if sym != inv(sym))
            symInv = inv(sym);
            d = angle(mori, symInv);
            if d < minAng
                minAng = d;
                bestRot = symInv;
            end
        end
        
        % Store rotation for this grain
        grainRotations(bid) = bestRot;
    end
    
    % Apply rotations to EBSD data in one go
    validPixel = ebsd.grainId > 0;
    pixelRotations = grainRotations(ebsd.grainId(validPixel));
    toUpdate = ~isnan(pixelRotations);
    
    if any(toUpdate)
        % Get indices into ebsd that need updating
        validIndices = find(validPixel);
        updateIndices = validIndices(toUpdate);
        
        ebsd(updateIndices).orientations = ebsd(updateIndices).orientations .* pixelRotations(toUpdate);
        countChanged = length(updateIndices);
    end
    
    fprintf('Updated %d EBSD data points (from %d merged grains).\n', countChanged, length(uBadIds));

    % 6. Merge the grains
    % This updates the grain structure
    [grains, parentId] = merge(grains, gB_to_merge);
    
    % Update EBSD grain IDs to match the merged result
    ebsd.grainId = parentId(ebsd.grainId);
    
end
