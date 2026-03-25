%% batch_GBCD_Analysis.m
% Batch script to calculate GBCD/GBND from multiple .h5oina datasets, save the 
% results individually, and calculate the weighted average of all datasets.

clear; close all; clc;
plotx2west;
plotzOutOfPlane;

%% 1. Setup Directories and Options
do_GBCD = false; % Calculate GBCD for specific sigma boundaries
do_GBND = true; % Calculate GBND for all boundaries

dataDir = 'R:\Scratch\201\Annalena Erlacher\3D_Mg_1750_4';
files = dir(fullfile(dataDir, '*.h5oina'));

if isempty(files)
    error('No .h5oina files found in the specified directory.');
end

% Allow user to interactively select/deselect files
fileNames = {files.name};
[indx, tf] = listdlg('PromptString', 'Select files to process:', ...
                     'SelectionMode', 'multiple', ...
                     'ListString', fileNames, ...
                     'InitialValue', 1:length(fileNames), ...
                     'ListSize', [300, 300]);

if tf == 0
    disp('No files selected. Exiting script.');
    return;
end

files = files(indx);
numFiles = length(files);
fprintf('Found %d datasets to process.\n', numFiles);

%% 2. Loop through datasets and calculate Properties
if do_GBCD
    gbcd_list = cell(numFiles, 1);
    gbcd_len_list = zeros(numFiles, 1);
end
if do_GBND
    gbnd_list = cell(numFiles, 1);
    gbnd_len_list = zeros(numFiles, 1);
end
CS_ref = [];

for i = 1:numFiles
    fname = files(i).name;
    fullPath = fullfile(dataDir, fname);
    [~, baseName, ~] = fileparts(fname);
    matFile = fullfile(dataDir, [baseName, '_GBCD.mat']);
    matFile_GBND = fullfile(dataDir, [baseName, '_GBND.mat']);
    
    fprintf('\n--- Processing File %d/%d: %s ---\n', i, numFiles, fname);
    
    need_calc_GBCD = do_GBCD && ~exist(matFile, 'file');
    need_calc_GBND = do_GBND && ~exist(matFile_GBND, 'file');
    
    % Check if GBCD .mat already exists
    if do_GBCD && ~need_calc_GBCD
        fprintf('  -> Loading existing GBCD results from %s\n', matFile);
        loadData = load(matFile, 'gbcd', 'len', 'CS');
        gbcd_list{i} = loadData.gbcd;
        gbcd_len_list(i) = loadData.len;
        if isempty(CS_ref)
            CS_ref = loadData.CS;
        end
    end
    
    % Check if GBND .mat already exists
    if do_GBND && ~need_calc_GBND
        fprintf('  -> Loading existing GBND results from %s\n', matFile_GBND);
        loadData = load(matFile_GBND, 'gbnd', 'len', 'CS');
        gbnd_list{i} = loadData.gbnd;
        gbnd_len_list(i) = loadData.len;
        if isempty(CS_ref)
            CS_ref = loadData.CS;
        end
    end
    
    % Perform calculations if missing
    if need_calc_GBCD || need_calc_GBND
        fprintf('  -> Loading and reconstructing map data...\n');
        ebsd_raw = loadEBSD_h5oina(fullPath, 'convertEuler2SpatialReferenceFrame');
        CS = ebsd_raw.CS;
        if isempty(CS_ref), CS_ref = CS; end
        
        % Perform corrections for misaligned axes
        rot = rotation.byAxisAngle(yvector, 180*degree);
        ebsd_raw = rotate(ebsd_raw, rot, 'keepEuler');
        ebsd_raw.orientations = project2FundamentalRegion(ebsd_raw.orientations);
        
        % Reconstruct grains once locally to save time if both are selected
        pseudoSym1 = orientation.byAxisAngle(Miller(0,0,0,1,CS), 60*degree);
        pseudoSym2 = orientation.byAxisAngle(Miller(0,0,0,1,CS), 30*degree);
        pseudoSym = [pseudoSym1, pseudoSym2];
        
        ebsd_clean = pseudoSymmetryCorrection(ebsd_raw, pseudoSym);
        [grains, ebsd_clean.grainId] = calcGrains(ebsd_clean, 'angle', 5*degree, 'alpha', 1, 'minPixel', 5);
        grains = smooth(grains, 5);
        gB = grains.boundary('indexed', 'indexed');
        
        if need_calc_GBCD
            fprintf('  -> Calculating GBCD...\n');
            sigmaAngles = [60, 38.21, 27.8, 46.83, 21.79] * degree;
            sigmas = orientation.byAxisAngle(Miller(0,0,0,1,CS), sigmaAngles);
            [gbcd, len] = performGBCDAnalysis(gB, grains, sigmas);
            
            gbcd_list{i} = gbcd;
            gbcd_len_list(i) = len;
            save(matFile, 'gbcd', 'len', 'CS');
        end
        
        if need_calc_GBND
            fprintf('  -> Calculating Total GBND...\n');
            [gbnd, len] = performGBCDAnalysis(gB, grains);
            gbnd_list{i} = gbnd;
            gbnd_len_list(i) = len;
            save(matFile_GBND, 'gbnd', 'len', 'CS');
        end
    end
end

%% 3. Calculate Weighted Averages & Plot Results
if do_GBCD
    fprintf('\nCalculating weighted average of GBCD datasets...\n');
    total_length_gbcd = sum(gbcd_len_list);
    gbcd_avg = gbcd_list{1} * gbcd_len_list(1);
    for i = 2:numFiles
        gbcd_avg = gbcd_avg + (gbcd_list{i} * gbcd_len_list(i));
    end
    gbcd_avg = gbcd_avg / total_length_gbcd;
    
    figure; plot(gbcd_avg);
    title(sprintf('Weighted Average GBCD (%d Datasets)', numFiles));
    mtexColorbar;
    hold on;
    annotate(Miller(0,0,0,1, CS_ref), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    annotate(Miller(1,1,-2,0, CS_ref), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    annotate(Miller(-1,1,0,0, CS_ref), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    hold off; drawNow(gcm, 'figSize', 'large');
    
    save(fullfile(dataDir, 'Weighted_Average_GBCD.mat'), 'gbcd_avg', 'total_length_gbcd', 'CS_ref');
end

if do_GBND
    fprintf('\nCalculating weighted average of GBND datasets...\n');
    total_length_gbnd = sum(gbnd_len_list);
    gbnd_avg = gbnd_list{1} * gbnd_len_list(1);
    for i = 2:numFiles
        gbnd_avg = gbnd_avg + (gbnd_list{i} * gbnd_len_list(i));
    end
    gbnd_avg = gbnd_avg / total_length_gbnd;
    
    figure; plot(gbnd_avg);
    title(sprintf('Weighted Average GBND (%d Datasets)', numFiles));
    mtexColorbar;
    hold on;
    annotate(Miller(0,0,0,1, CS_ref), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    annotate(Miller(1,1,-2,0, CS_ref), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    annotate(Miller(-1,1,0,0, CS_ref), 'labeled', 'backgroundColor', 'w', 'fontWeight', 'bold');
    hold off; drawNow(gcm, 'figSize', 'large');
    
    save(fullfile(dataDir, 'Weighted_Average_GBND.mat'), 'gbnd_avg', 'total_length_gbnd', 'CS_ref');
end

fprintf('\nFinished processing!\n');