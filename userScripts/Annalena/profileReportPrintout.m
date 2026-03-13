% HOW TO USE:
% 1. Run `profile on` in the command window.
% 2. Execute your MTEX script as you normally would.
% 3. Run `profile off`.
% 4. Update the `funcName` variable below and run this block of code.
% 5. Optional: You can also view the results by running: profile viewer

funcName = 'graphedPseudoSymRemoval'; % Set your target function here

p = profile('info');
idx = endsWith({p.FunctionTable.FunctionName}, funcName);

if ~any(idx)
    error('Function "%s" not found in profiler data. Did it run?', funcName);
end

fData = p.FunctionTable(idx);
fData = fData(1); % Take the first match if there are multiple

% Find the directory of the function to save the MD file there
[funcPath, ~, ~] = fileparts(which(funcName));
if isempty(funcPath)
    funcPath = pwd; % Fallback to current working directory
end

% Create dynamic filename with date
dateStr = char(datetime('now', 'Format', 'yyyy-MM-dd_hh-mm-ss'));
fileName = fullfile(funcPath, sprintf('%s_profiler_%s.md', funcName, dateStr));

% Write the Markdown file
fid = fopen(fileName, 'w');
fprintf(fid, '# Profiler Report: %s\n', fData.FunctionName);
fprintf(fid, '**Total Time:** %.3fs | **Total Calls:** %d\n\n', fData.TotalTime, fData.NumCalls);
fprintf(fid, '| Line Number | Time (s) | Number of Calls |\n');
fprintf(fid, '| :--- | :--- | :--- |\n');

for i = 1:size(fData.ExecutedLines, 1)
    fprintf(fid, '| %d | %.4f | %d |\n', ...
        fData.ExecutedLines(i,1), fData.ExecutedLines(i,3), fData.ExecutedLines(i,2));
end

fclose(fid);
fprintf('Successfully exported profiler report to:\n%s\n', fileName);