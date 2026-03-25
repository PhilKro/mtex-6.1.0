function exportScaledFigure(fig, fname, varargin)
% exportScaledFigure Export a figure with scaled shorter side
%
%   exportScaledFigure(fig, fname)
%       exports the figure FIG to file FNAME with default settings:
%       (shorter side = 1300 px at 300 dpi), and appends a datetime string
%       before the file extension.
%
%   exportScaledFigure(fig, fname, 'target_px', N, 'dpi', D, 'noDateString')
%       allows you to override target size and dpi, and disable the date string.
%
%   Example:
%       f = figure;
%       plot(rand(10,1));
%       exportScaledFigure(f, 'myplot.png') 
%       % → myplot_09_04_15_30_12.png
%
%       exportScaledFigure(f, 'myplot.png', 'noDateString')
%       % → myplot.png

    % --- Defaults ---
    target_px = 1500;
    dpi = 300;
    addDateString = true;
    extraheight = false;
    extrawidth = false;
    % --- Parse optional inputs ---
    k = 1;
    while k <= numel(varargin)
        switch lower(varargin{k})
            case 'target_px'
                target_px = varargin{k+1};
                k = k + 1;
            case 'dpi'
                dpi = varargin{k+1};
                k = k + 1;
            case 'nodatestring'
                addDateString = false;
            case 'extraheight'
                extraheight = true;
            case 'extrawidth'
                extrawidth = true;
            otherwise
                error('Unknown option: %s', varargin{k});
        end
        k = k + 1;
    end

    % --- Add datetime string if enabled ---
    [filepath, name, ext] = fileparts(fname);
    if isempty(ext)
        ext = '.png'; % default extension
    end
    if addDateString
        dtstr = string(datetime('now','Format','MM_dd_HH_mm_ss'));
        name = name + "_" + dtstr;
    end
    fname = fullfile(filepath, append(name, ext));

    % --- Compute target size in inches ---
    target_in = target_px / dpi;

    % --- Read current figure size in inches ---
    oldUnits = fig.Units;
    fig.Units = 'inches';
    pos = fig.Position;   % [x, y, width, height] in inches

    % --- Scale so shorter side = target_in ---
    shorter_side = min(pos(3:4));
    scale = target_in / shorter_side;
    if extraheight
        pos(4)= pos(4)*1.3;
    elseif extrawidth
        pos(3)= pos(3)*1.3;
    end
    fig.Position(3:4) = pos(3:4) * scale;

    % --- Restore units ---
    fig.Units = oldUnits;

    % --- Export ---
    exportgraphics(fig, fname, 'Resolution', dpi);

    % --- Report actual export size ---
    width_px  = round(fig.Position(3) * dpi);
    height_px = round(fig.Position(4) * dpi);
    fprintf('Exported %s at %dx%d px (%d dpi)\n', fname, width_px, height_px, dpi);
end
