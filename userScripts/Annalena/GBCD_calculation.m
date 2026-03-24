%% Initial Setup
clear
close all
plotx2west
plotzOutOfPlane

%% 1. Load Data

% crystal symmetry, make sure to use correct space group and c/a ratio, currently Zn is
% used
dir = 'userScripts/Annalena/';
fname = 'AA-04 500 Nd 1500 Specimen 1 Site 5 Map Data 4.h5oina';
path = append(dir, fname);

%if loading the CS from the h5oina file
ebsd_raw = loadEBSD_h5oina(path, 'convertEuler2SpatialReferenceFrame');
CS = ebsd_raw.CS;

%smallest possible testSet
% ebsd_raw = ebsd_raw(inpolygon(ebsd_raw, [5 2 5 5]));

%small testSet includes straight boundary
% ebsd_raw = ebsd_raw(inpolygon(ebsd_raw, [5 2 50 30]));

%larger testSet includes this coral grains
% ebsd_raw = ebsd_raw(inpolygon(ebsd_raw, [5 2 50 200]));

%perform corrections for misaligned axes
rot = rotation.byAxisAngle(yvector,180*degree);
ebsd_raw = rotate(ebsd_raw,rot,'keepEuler');
ebsd_raw.orientations = project2FundamentalRegion(ebsd_raw.orientations);

%% test alignment of crystal reference frame with specimen coordinates
% use this colorkey for comparison with oxford
% ipf=ipfTSLKey(ebsd.CS);
% this is a better suited colorkey
ipf = ipfColorKey(ebsd_raw.CS);

% figure;
% plot(ipf)
% figure;
% plot(ebsd.CS,'figSize','small')
% annotate(ebsd.CS.aAxis,'MarkerFaceColor','r','label','a','backgroundColor','w')
% annotate(ebsd.CS.bAxis,'MarkerFaceColor','r','label','b','backgroundColor','w')
% annotate(-vector3d.Y,'MarkerFaceColor','green','label','-y','backgroundColor','w')
% annotate(vector3d.X,'MarkerFaceColor','green','label','x','backgroundColor','w')

%% plot ebsd
% figure;
% ipf.inversePoleFigureDirection=xvector;
% colors = ipf.orientation2color(ebsd('indexed').orientations);
% plot(ebsd('indexed'),colors)
% title('IPF X')
% 
% ipf.inversePoleFigureDirection = yvector;
% colors = ipf.orientation2color(ebsd('indexed').orientations);
% nextAxis
% plot(ebsd('indexed'),colors)
% title('IPF Y')
% exportScaledFigure(gcf, 'IPF_X_Y.png')

% ipf.inversePoleFigureDirection = zvector;
% colors = ipf.orientation2color(ebsd('indexed').orientations);
% nextAxis
% plot(ebsd('indexed'),colors)
% title('IPF Z')

%% Initial grain reconstruction and pseudo symmetry cleanup
% 1. Initial rough grain calc

[grains_raw, ebsd_raw.grainId] = calcGrains(ebsd_raw, 'angle', 5*degree, 'alpha',1, 'minPixel',5);

% figure; plot(grains, grains.meanOrientation, 'noBoundary'); title('Pre-cleaned grains');
% % exportScaledFigure(gcf, 'Pre-cleaned grains.jpg', 'dpi', 200)
%% plot sigma boundaries befor cleanup
plotgrains = grains_raw;


sigmaAngles = [60, 38.21, 27.8, 46.83, 21.79] * degree;
sigmaNames = {'Sigma 3', 'Sigma 7', 'Sigma 13', 'Sigma 19', 'Sigma 21'};
colors = {'r', 'g', 'b', 'c', 'm'};
grainColors = ipf.orientation2color(plotgrains.meanOrientation);
plotgBs = plotgrains.boundary('indexed','indexed');

sigmas = orientation.byAxisAngle(Miller(0,0,0,1,CS), sigmaAngles);

wierd_sigmas = orientation.byAxisAngle(Miller(1,-2,1,0,CS,'UVWT'), 20*degree);

isSigmaAny = false(size(gB_raw));

% Visualize where they are
figure;
plot(plotgrains, grainColors, 'faceAlpha', 0.1, 'noBoundary')
hold on

for i = 1:length(sigmas)
    isSigma = angle(plotgBs.misorientation, sigmas(i)) < 5*degree;
    isSigmaAny = isSigmaAny | isSigma;
    
    plot(plotgBs(isSigma), 'lineColor', colors{i}, 'lineWidth', 2, 'DisplayName', sigmaNames{i})
end

legend show
title('Identified Sigma Boundaries')
hold off

clear plotgrains
clear grainColors
clear plotgBs

%% 3. Define the pseudo symmetry (180 deg rotation around c-axis)
pseudoSym1 = orientation.byAxisAngle(Miller(0,0,0,1,CS), 60*degree);
pseudoSym2 = orientation.byAxisAngle(Miller(0,0,0,1,CS), 30*degree);
pseudoSym = [pseudoSym1, pseudoSym2];

%% 3. PseudoSymmetry removal

% Version 1: Ralph orginal (dauert zuuu lange)
% [ebsd,grains] = cleanUpPseudoSym(ebsd_raw,grains_raw,pseudoSym(1));

% Version 2: Ralph aber schnell (~6000 pixel verändert)
% veränder den threshold parameter
% [ebsd,grains, numChanged] = cleanUpPseudoSym_Phil(ebsd_raw,grains_raw,pseudoSym, 'threshold', 1.5);
% numChanged_phi

% Version 3: Philipp (26000 pixel verändert)
% ebsd = pseudoSymmetryCorrection(ebsd_raw, pseudoSym);

%% Parameter Analysis for Pseudo Symmetry Grains (Debug)
% % Identify grains involved in pseudo-symmetry boundaries
% %grains = smooth(grains, 1); 
% 
% gB_check = grains.boundary('indexed','indexed');
% isPseudo_check = false(size(gB_check));
% for k = 1:length(pseudoSym)
%     isPseudo_check = isPseudo_check | (angle(gB_check.misorientation, pseudoSym(k)) < 5*degree);
% end
% 
% grainIds_check = unique(gB_check(isPseudo_check).grainId);
% grainIds_check(grainIds_check == 0) = [];
% grains_check = grains(grainIds_check);
% 
% check_curv = nan(length(grains_check),1);
% check_PA = nan(length(grains_check),1);
% check_sz = nan(length(grains_check),1);
% 
% for i = 1:length(grains_check)
%    check_curv(i) = mean(abs(grains_check(i).boundary.curvature));
%    check_PA(i) = sum(grains_check(i).boundary.segLength) / sum(grains_check(i).numPixel);
%    check_sz(i) = sum(grains_check(i).numPixel);
% end
% 
% figure; scatter(check_PA, check_curv); xlabel('sum(segLength)/numPixel'); ylabel('Curvature'); title('Curvature vs Surface/Volume');
% figure; scatter(check_sz, check_PA); xlabel('numPixel'); ylabel('sum(segLength)/numPixel'); title('Surface/Volume vs Size');
% figure; plot(grains.boundary); hold on; plot(grains_check, check_curv); mtexColorbar
% figure; plot(grains.boundary); hold on; plot(grains_check, check_PA); mtexColorbar
%% Final grain reconstruction

% [grains, ebsd('indexed').grainId] = calcGrains(ebsd('indexed'), 'angle', 5*degree, 'alpha',1, 'minPixel',5);

% Smooth boundaries to get better trace directions
grains = smooth(grains, 5);
grainColors = ipf.orientation2color(grains.meanOrientation);

% Extract all grain boundaries
gB = grains.boundary('indexed','indexed');

figure; plot(grains, grains.meanOrientation, 'noBoundary'); title('Full clean')
% exportScaledFigure(gcf, 'CleanedPseudoSym_Grains.jpg', 'dpi', 200)

%% 2. Plot sigma boundaries after cleanup
plotgrains = grains;


sigmaAngles = [60, 38.21, 27.8, 46.83, 21.79] * degree;
sigmaNames = {'Sigma 3', 'Sigma 7', 'Sigma 13', 'Sigma 19', 'Sigma 21'};
colors = {'r', 'g', 'b', 'c', 'm'};
grainColors = ipf.orientation2color(plotgrains.meanOrientation);
plotgBs = plotgrains.boundary('indexed','indexed');

sigmas = orientation.byAxisAngle(Miller(0,0,0,1,CS), sigmaAngles);

wierd_sigmas = orientation.byAxisAngle(Miller(1,-2,1,0,CS,'UVWT'), 20*degree);

isSigmaAny = false(size(gB_raw));

% Visualize where they are
figure;
plot(plotgrains, grainColors, 'faceAlpha', 0.1, 'noBoundary')
hold on

for i = 1:length(sigmas)
    isSigma = angle(plotgBs.misorientation, sigmas(i)) < 5*degree;
    isSigmaAny = isSigmaAny | isSigma;
    
    plot(plotgBs(isSigma), 'lineColor', colors{i}, 'lineWidth', 2, 'DisplayName', sigmaNames{i})
end

legend show
title('Identified Sigma Boundaries')
hold off

clear plotgrains
clear grainColors
clear plotgBs

%% Plot grains with crystal shapes
% cS = crystalShape.hex(ebsd.CS);
% 
% figure;
% ipf.inversePoleFigureDirection = zvector;
% grainColors = ipf.orientation2color(grains.meanOrientation);
% 
% plot(grains,grainColors,'FaceAlpha',0.5,'linewidth',2)
% isBig = grains.numPixel>100;
% 
% % define a list of crystal shape that is oriented as the grain mean
% % orientation and scaled according to the grain area
% cSGrains = grains(isBig).meanOrientation * cS * 0.7 * sqrt(grains(isBig).area);
% 
% % now we can plot these crystal shapes at the grain centers
% hold on
% plot(grains(isBig).centroid + cSGrains,'FaceColor','r','FaceAlpha',0.7)
% hold off
% drawNow(gcm,'final')
% exportScaledFigure(gcf, 'Grains_with_CrystalShapes.png')

%% 3. Identify Sigma 3 Boundaries (The 3 Misorientation Parameters)

% Define Sigma misorientations (axis <0001>)


isSigmaAny = false(size(gB));

% Visualize where they are
figure;
plot(grains, grainColors, 'faceAlpha', 0.1, 'noBoundary')
hold on

for i = 1:length(sigmas)
    isSigma = angle(gB.misorientation, sigmas(i)) < 5*degree;
    isSigmaAny = isSigmaAny | isSigma;
    
    plot(gB(isSigma), 'lineColor', colors{i}, 'lineWidth', 2, 'DisplayName', sigmaNames{i})
end

legend show
title('Identified Sigma Boundaries')
hold off
% exportScaledFigure(gcf, 'Identified_Sigma_Boundaries.jpg', 'dpi', 200)

%% 4. Calculate Grain Boundary Normal Distribution

% GBND for ALL boundaries
gbnd = calcGBND(gB,grains,'halfwidth',7.5*degree);

% GBCD for sigma boundaries
gbcd = calcGBND(gB(isSigmaAny), grains,'halfwidth',7.5*degree);

% Loop through the sigmas and add up the spherical harmonics afterwards
gbcd_reccomended = 0;
for i = 1:length(sigmas)
    gbcd_reccomended = gbcd_reccomended + calcGBND(gB,grains,sigmas(i),'halfwidth',7.5*degree);
end

%% 5. Visualize the 5-Parameter Character
figure;

% Plot the GBND in the fundamental region of the crystal
plot(gbnd) %, 'fundamentalRegion')
mtexColorbar
title('GBND (Crystal Frame)')

% Annotate the expected boundary plane {0001}
hold on
h = Miller(0,0,0,1, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
h = Miller(1,1,-2,0, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
h = Miller(-1,1,0,0, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
% exportScaledFigure(gcf, 'GBND_CrystalFrame.jpg')

figure;

% Plot the GBCD in the fundamental region of the crystal
plot(gbcd) %, 'fundamentalRegion')
mtexColorbar
title('GBCD for sigma boundaries (Crystal Frame)')

% Annotate the expected boundary plane {0001}
hold on
h = Miller(0,0,0,1, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
h = Miller(1,1,-2,0, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
h = Miller(-1,1,0,0, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
% exportScaledFigure(gcf, 'GBCD_Sigma_CrystalFrame.jpg')

figure;

% Plot the GBCD, calculated the recommended way
plot(gbcd) %, 'fundamentalRegion')
mtexColorbar
title('GBCD (recommended) for sigma boundaries (Crystal Frame)')

% Annotate the expected boundary plane {0001}

hold on
h = Miller(0,0,0,1, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
h = Miller(1,1,-2,0, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
h = Miller(-1,1,0,0, CS);
annotate(h, 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
% exportScaledFigure(gcf, 'GBCD_Recommended_Sigma_CrystalFrame.jpg')
%% 6. Checking consistency for individual segments
% While we can't know the exact plane for a single segment, we can check
% if the observed trace is CONSISTENT with a specific plane .

sigmaGBS = gB(isSigmaAny);

[~,gbcd_max_miller] = max(gbcd);

traces = sigmaGBS.direction;

grainsnext2traces = grains(sigmaGBS.grainId(:,2));

traces_in_GrainCoords = inv(grainsnext2traces.meanOrientation).*traces;

angletrace2maxgbcd = angle(traces_in_GrainCoords, gbcd_max_miller)/degree;
% hist(angletrace2maxgbcd,100)

%fractionCoherent = sum(gB3.segLength(isCoherent)) / sum(gB3.segLength);

%fprintf('Fraction of Sigma 3 boundaries consistent with coherent {111} plane: %.1f%%\n', fractionCoherent * 100);

% Highlight coherent vs incoherent candidates
figure;
plot(grains, grainColors, 'faceAlpha', 0.15, 'noBoundary')
hold on
quiver(grains(sigmaGBS.grainId(:,1)), grainsnext2traces.meanOrientation.*gbcd_max_miller, 'color','r') % cross(grainsnext2traces.meanOrientation.*gbcd_max_miller, zvector)
plot(sigmaGBS, angletrace2maxgbcd, 'lineWidth', 2, 'DisplayName', 'Inconsistent with GBCD max')
mtexColorbar
%plot(sigmaGBS(trace_in_maxgbcd), 'lineColor', 'r', 'lineWidth', 2, 'DisplayName', 'Consistent with GBCD max')
legend
title('Trace Consistency Analysis')
% exportScaledFigure(gcf, 'Trace_Consistency_Analysis.jpg', 'dpi', 200)

%% GBCD averaging from different maps

current_gbcd = gbcd;
current_boundary_length = sum(gB(isSigmaAny).segLength);

% 1. Display the required warning to the user
warning(['Only one GBCD is calculated in this script. ', ...
         'To calculate a weighted average, the second GBCD must be imported ', ...
         'from a .mat file. Please ensure the .mat file contains variables ', ...
         'appropriately named ''imported_gbcd'' and ''imported_length''.']);

% 2. Prompt the user to select and load the .mat file
[fileName, pathName] = uigetfile('*.mat', 'Select the .mat file containing the second GBCD');

if isequal(fileName, 0)
    disp('File selection canceled. Skipping the weighted average GBCD calculation.');
else
    % Load the data from the selected .mat file into a struct
    loadedData = load(fullfile(pathName, fileName));
    
    % Verify the appropriately named variables exist in the loaded data
    if isfield(loadedData, 'imported_gbcd') && isfield(loadedData, 'imported_length')
        
        imported_gbcd = loadedData.imported_gbcd;
        imported_length = loadedData.imported_length;
        
        % Ensure the current and imported GBCDs have identical matrix/array sizes
        if isequal(size(current_gbcd), size(imported_gbcd))
            
            % 3. Calculate the weighted average based on boundary lengths
            total_length = current_boundary_length + imported_length;
            
            weighted_avg_gbcd = ((current_gbcd * current_boundary_length) + ...
                                 (imported_gbcd * imported_length)) / total_length;
                                 
            disp('Success: Weighted average GBCD calculated and saved to variable ''weighted_avg_gbcd''.');
            
        else
            error('Dimension mismatch: The imported GBCD array does not match the size of the calculated GBCD array.');
        end
    else
        error('Missing variables: The selected .mat file must contain exactly ''imported_gbcd'' and ''imported_length''.');
    end
end
