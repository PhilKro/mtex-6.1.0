function [gbcd, total_len] = performGBCDAnalysis(data, ref_data, varargin)
    % PERFORMGBCDANALYSIS calculates the Grain Boundary Character/Normal Distribution
    %
    % Usage 1 (from pre-calculated grains and boundaries):
    %   [gbcd, len] = performGBCDAnalysis(gB, grains, [sigmas])

    sigmas = getClass(varargin, 'orientation');


    if isa(data, 'grainBoundary')
        gB = data;
        grains = ref_data;
    else
        error('First argument must be of type grainBoundary');
    end
    
    % 2. Check for optional sigma input
    if ~isempty(sigmas)
        % Evaluate sigma boundary properties for lengths
        isSigmaAny = false(size(gB));
        for i = 1:length(sigmas)
            isSigmaAny = isSigmaAny | (angle(gB.misorientation, sigmas(i)) < 5*degree);
        end
        target_gB = gB(isSigmaAny);
    else
        % No sigma input, use all boundaries
        target_gB = gB;
    end
    
    total_len = sum(target_gB.segLength);
    
    % 3. Compute GBCD/GBND
    gbcd = calcGBND(target_gB, grains, 'halfwidth', 7.5*degree);
end