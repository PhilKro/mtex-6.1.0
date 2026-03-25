function [gbcd, total_len] = performGBCDAnalysis(data, ref_data, varargin)
    % PERFORMGBCDANALYSIS calculates the Grain Boundary Character/Normal Distribution
    %
    % Usage 1 (from raw EBSD):
    %   [gbcd, len] = performGBCDAnalysis(ebsd_data, CS, [sigmas])
    %
    % Usage 2 (from pre-calculated grains and boundaries):
    %   [gbcd, len] = performGBCDAnalysis(gB, grains, [sigmas])

    sigmas = getClass(varargin, 'orientation');

    if isa(data, 'EBSD')
        ebsd_data = data;
        CS = ref_data;
        
        % 1. Reconstruct grains and cleanup pseudo-symmetry
        pseudoSym1 = orientation.byAxisAngle(Miller(0,0,0,1,CS), 60*degree);
        pseudoSym2 = orientation.byAxisAngle(Miller(0,0,0,1,CS), 30*degree);
        pseudoSym = [pseudoSym1, pseudoSym2];
        
        % Initial grain reconstruction
        [grains_raw, ebsd_data.grainId] = calcGrains(ebsd_data, 'angle', 5*degree, 'alpha', 1, 'minPixel', 5);
        
        % Cleanup with cleanUpPseudoSym_Phil
        [ebsd_data, grains, numChanged] = cleanUpPseudoSym_Phil(ebsd_data, grains_raw, pseudoSym, 'threshold', 1.5);
        disp(['Number of pixels changed during pseudo-symmetry cleanup: ', num2str(numChanged)]);
        grains = smooth(grains, 5);
        gB = grains.boundary('indexed', 'indexed');
    elseif isa(data, 'grainBoundary')
        gB = data;
        grains = ref_data;
    else
        error('First argument must be of type EBSD or grainBoundary');
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