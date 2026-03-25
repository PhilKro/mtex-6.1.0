% wenn du subsets von deinem ebsd datensatz hersellen willst (plottet schneller zum testen von reconstructionseinstellungen)
rect = [xstart, ystart, dx, dy];
ebsd = ebsd(inpolygon(ebsd_raw, rect));

% um ebsd daten auszuschliessen mit hohen MAD werten
ebsd(ebsd.prop.mad > 0.1) = 'notIndexed';

% gridification hilft um die daten schneller zu plotten
ebsd=ebsd.gridify;


% find other parameters type "edit calcGrains" to open the function. https://mtex-toolbox.github.io/EBSD.calcGrains.html, https://mtex-toolbox.github.io/GrainReconstruction.html
%%% perform reconstruction %%%
ebsd = ebsd0;

grainThresh = 5*degree; 
% min grainsize, in pixels
minSize = 3;
% aplha parameters, defines how much to fill during recons
% truction
alpha = 1;

% wenn du mehr lücken füllen willst, nur die indizierten ebsd daten verwenden
[grains,ebsd('indexed').grainId,ebsd('indexed').mis2mean] = calcGrains(ebsd('indexed'),'minPixel',minSize,'angle',grainThresh,'alpha',alpha);


% klassische grain reconstruction
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);

% --- machnmal hilfreich, wenn zu viele kleine Körner da sind (für dich warscheinlich nicht relecant tatsächlich)
ebsd(grains(grains.numPixel<=minSize)) = 'notIndexed';
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'minPixel',minSize,'angle',grainThresh,'alpha',alpha);
% --- 