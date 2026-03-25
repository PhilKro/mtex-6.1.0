%% mapSegmentation_test.m
% Takes the EBSD map of the alumina and runs the GBCD analysis (with cleanup) 
% on the map. Then takes the same map, segments it into 4 pieces using inpolygon,
% runs the GBCD analysis (with cleanup) on the 4 pieces, calculates the weighted 
% (by segment length) average of the 4 GBCDs, and compares it to the GBCD of the whole map.

clear
close all
plotx2west
plotzOutOfPlane

%% 1. Load Data
dir = 'C:\Users\phkr\OneDrive - empa.ch\mtex-6.1.0\userScripts\Annalena';
fname = 'AA-04 500 Nd 1500 Specimen 1 Site 5 Map Data 4.h5oina';
path = fullfile(dir, fname);

ebsd_raw = loadEBSD_h5oina(path, 'convertEuler2SpatialReferenceFrame');
CS = ebsd_raw.CS;

% ebsd_raw = ebsd_raw(inpolygon(ebsd_raw, [0 0 50 30]));

% Perform corrections for misaligned axes
rot = rotation.byAxisAngle(yvector, 180*degree);
ebsd_raw = rotate(ebsd_raw, rot, 'keepEuler');
ebsd_raw.orientations = project2FundamentalRegion(ebsd_raw.orientations);

sigmaAngles = [60, 38.21, 27.8, 46.83, 21.79] * degree;
sigmas = orientation.byAxisAngle(Miller(0,0,0,1,CS), sigmaAngles);

%% 2. GBCD Analysis on the WHOLE Map
fprintf('\n--- Processing WHOLE map ---\n');
[gbcd_whole, len_whole] = performGBCDAnalysis(ebsd_raw, CS, sigmas);

%% 3. Segment the Map into 4 Pieces using inpolygon
xmin = min(ebsd_raw.x); xmax = max(ebsd_raw.x);
ymin = min(ebsd_raw.y); ymax = max(ebsd_raw.y);
xmid = (xmin + xmax) / 2;
ymid = (ymin + ymax) / 2;

% Calculate extents for the quadrants
dx = xmid - xmin;
dy = ymid - ymin;

% Define rectangles for the 4 quadrants using [xmin, ymin, dx, dy] syntax
% "-" are needed because of the 'keepEuler' rotation, much cursed
rect1 = [xmin, ymin, -dx, dy];
rect2 = [-xmid, ymin, -dx, dy];
rect3 = [xmin, ymid, -dx, dy];
rect4 = [-xmid, ymid, -dx, dy];

ebsd_q1 = ebsd_raw(inpolygon(ebsd_raw, rect1));
ebsd_q2 = ebsd_raw(inpolygon(ebsd_raw, rect2));
ebsd_q3 = ebsd_raw(inpolygon(ebsd_raw, rect3));
ebsd_q4 = ebsd_raw(inpolygon(ebsd_raw, rect4));

%% 4. GBCD Analysis on the 4 Pieces
fprintf('\n--- Processing Quadrant 1 ---\n');
[gbcd_q1, len_q1] = performGBCDAnalysis(ebsd_q1, CS, sigmas);
% save("%s.mat", dataensatznummer, "gbcd_q1", "len_q1"); % Debug save disabled

fprintf('\n--- Processing Quadrant 2 ---\n');
[gbcd_q2, len_q2] = performGBCDAnalysis(ebsd_q2, CS, sigmas);

fprintf('\n--- Processing Quadrant 3 ---\n');
[gbcd_q3, len_q3] = performGBCDAnalysis(ebsd_q3, CS, sigmas);

fprintf('\n--- Processing Quadrant 4 ---\n');
[gbcd_q4, len_q4] = performGBCDAnalysis(ebsd_q4, CS, sigmas);

%% 5. Calculate Weighted Average
fprintf('\nCalculating weighted average of the 4 quadrants...\n');
total_segment_length = len_q1 + len_q2 + len_q3 + len_q4;

gbcd_avg = (gbcd_q1 * len_q1 + ...
            gbcd_q2 * len_q2 + ...
            gbcd_q3 * len_q3 + ...
            gbcd_q4 * len_q4) / total_segment_length;

%% 6. Plot and Compare Results
figure;
newMtexFigure('layout', [1, 2]);

% Plot 1: Whole Map GBCD
nextAxis(1, 1);
plot(gbcd_whole);
title('GBCD (Whole Map)');
annotateExpectedPlanes(CS);

% Plot 2: Weighted Average GBCD
nextAxis(1, 2);
plot(gbcd_avg);
title('GBCD (Weighted Avg of 4 Quadrants)');
annotateExpectedPlanes(CS);
mtexColorbar;

drawNow(gcm, 'figSize', 'large');

%% Helper Functions
function annotateExpectedPlanes(CS)
    hold on;
    annotate(Miller(0,0,0,1, CS), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    annotate(Miller(1,1,-2,0, CS), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    annotate(Miller(-1,1,0,0, CS), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    hold off;
end
