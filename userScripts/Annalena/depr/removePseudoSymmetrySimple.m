function [grains, ebsd] = removePseudoSymmetrySimple(ebsd, grains, pseudoSym, threshold, varargin)
% REMOVEPSEUDOSYMMETRYSIMPLE merges pseudo-symmetry artifacts based on boundary perimeter/area ratio.
%
% Usage:
%   [grains, ebsd] = removePseudoSymmetrySimple(ebsd, grains, pseudoSym, threshold)
%
% Inputs:
%   ebsd      - EBSD variable
%   grains    - grain variable
%   pseudoSym - list of pseudo-symmetry misorientations
%   threshold - (Optional) threshold for sum(segLength)/numPixel (default 0.38)
%
% Options:
%   'ignoreMAD' - bool, ignore MAD for direction decision (default false)

    if nargin < 4 || isempty(threshold)
        threshold = 0.38;
    end
    
    ignoreMAD = false;
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'ignoreMAD')
            ignoreMAD = varargin{i+1};
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

    % 2. Prepare Analysis
    % Get IDs of grains on both sides of the boundary
    grainIds = gB_ps.grainId;
    
    % Filter out outer boundaries (ID=0)
    validMask = all(grainIds > 0, 2);
    gB_ps = gB_ps(validMask);
    grainIds = grainIds(validMask, :);
    
    % Sort grain IDs to ensure consistent pairs
    grainIds = sort(grainIds, 2);
    
    % Unique pairs involved in pseudo-symmetry
    [uPairs, ~, ~] = unique(grainIds, 'rows');
    
    % Calculate MAD if available/requested (for deciding merge direction)
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
    
    % 3. Calculate Ratio and Decide Merges
    % Criterion: sum(segLength) / numPixel > threshold
    
    % Pre-calculate for involved grains to avoid re-calculation in loop
    uGrainsList = unique(uPairs(:));
    gRatio = nan(max(uGrainsList), 1);
    
    for i = 1:length(uGrainsList)
        gid = uGrainsList(i);
        if gid <= length(grains) && grains(gid).numPixel > 0
            % Ratio: Perimeter / Area (pixels)
            gRatio(gid) = sum(grains(gid).boundary.segLength) / grains(gid).numPixel;
        end
    end
    
    badIds = [];
    goodIds = [];
    
    for k = 1:size(uPairs, 1)
        id1 = uPairs(k, 1);
        id2 = uPairs(k, 2);
        
        r1 = gRatio(id1);
        r2 = gRatio(id2);
        
        is1Bad = r1 > threshold;
        is2Bad = r2 > threshold;
        
        if is1Bad || is2Bad
            % At least one needs removal. Decide direction.
            
            % Preference: Merge Bad -> Good.
            % If both Bad (or both Good, impossible here), merge "Worse" -> "Better".
            % "Worse" defined by MAD (if avail) or Size (smaller is worse).
            
            is1WorseQuality = false;
            
            if useMAD && ~isnan(grainMeanMAD(id1)) && ~isnan(grainMeanMAD(id2))
                if grainMeanMAD(id1) > grainMeanMAD(id2)
                    is1WorseQuality = true;
                elseif grainMeanMAD(id1) < grainMeanMAD(id2)
                    is1WorseQuality = false;
                else
                    % Tie-break size: smaller is worse
                    is1WorseQuality = grains(id1).numPixel < grains(id2).numPixel;
                end
            else
                % Size only: smaller is worse
                is1WorseQuality = grains(id1).numPixel < grains(id2).numPixel;
            end
            
            % Enforce ratio threshold:
            % If only one exceeds threshold, that one MUST be the 'bad' one.
            if is1Bad && ~is2Bad
                is1WorseQuality = true;
            elseif is2Bad && ~is1Bad
                is1WorseQuality = false;
            end
            % If both exceed, we stick to is1WorseQuality based on MAD/Size.
            
            if is1WorseQuality
                badIds = [badIds; id1];
                goodIds = [goodIds; id2];
            else
                badIds = [badIds; id2];
                goodIds = [goodIds; id1];
            end
        end
    end
    
    if isempty(badIds)
        return;
    end
    
    % 4. Select segments to merge
    pairsToMerge = sort([badIds, goodIds], 2);
    maskMerge = ismember(grainIds, pairsToMerge, 'rows');
    gB_to_merge = gB_ps(maskMerge);
    
    if isempty(gB_to_merge)
        return;
    end

    % 5. Rotate 'Bad' Grains before Merging
    % We rotate 'bad' grains to match 'good' neighbors orientation
    
    % Resolve conflicts: if a bad grain maps to multiple good ones, pick one.
    [uBadIds, idx] = unique(badIds);
    uGoodIds = goodIds(idx);
    
    grainRotations = orientation.nan(length(grains), 1, pseudoSym(1).CS, pseudoSym(1).CS);
    
    for k = 1:length(uBadIds)
        bid = uBadIds(k);
        gid = uGoodIds(k);
        
        oriBad = grains(bid).meanOrientation;
        oriGood = grains(gid).meanOrientation;
        
        mori = inv(oriBad) * oriGood;
        
        bestRot = orientation.nan;
        minAng = inf;
        
        for s = 1:length(pseudoSym)
            sym = pseudoSym(s);
            d = angle(mori, sym);
            if d < minAng
                minAng = d;
                bestRot = sym;
            end
            % Check inverse symmetry
            symInv = inv(sym);
            d = angle(mori, symInv);
            if d < minAng
                minAng = d;
                bestRot = symInv;
            end
        end
        grainRotations(bid) = bestRot;
    end
    
    % Apply rotations to EBSD
    validPixel = ebsd.grainId > 0;
    pixelRotations = grainRotations(ebsd.grainId(validPixel));
    toUpdate = ~isnan(pixelRotations);
    
    if any(toUpdate)
        validIndices = find(validPixel);
        updateIndices = validIndices(toUpdate);
        ebsd(updateIndices).orientations = ebsd(updateIndices).orientations .* pixelRotations(toUpdate);
        fprintf('Updated %d EBSD data points (from %d merged grains).\n', length(updateIndices), length(uBadIds));
    end

    % 6. Merge Grains
    [grains, parentId] = merge(grains, gB_to_merge);
    ebsd.grainId = parentId(ebsd.grainId);
    
end