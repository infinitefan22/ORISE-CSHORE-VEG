function res = organize_files(files, matchtext)
% this function will organize files based on shared titles 
% files = the entire list of files you want organized into separate groups
% based on title
% matchtext = the format you want the files organized (what in the shared
% title repeats (exact words or format) ) |   '^(HighDensity|LowDensity|Baseline)_[\w]+_[\w]+(_\w+)*(?=_Trial)'
fileGroups = struct();
for i = 1:length(files)
    filename = files(i).name; % Extract the filename
    
    % Use regular expression to extract the part before the trial number (e.g., HighDensity_Height1_Height2)
    parts = regexp(filename, matchtext, 'match');
    
    if ~isempty(parts)
        commonPart = parts{1} ; %making it just the string
        
        if ~any(strcmp(fieldnames(fileGroups), commonPart))  % Check if the common part already exists in the groups
        % If not found, create a new field for the category in the structure
            fileGroups.(commonPart) = {filename};  % Creates a field named "commonPart" and stores the filename in a cell array
        else
            % If found, add the file to the existing category field in the structure
            fileGroups.(commonPart){end+1} = filename;
        end
    end
end

res = fileGroups ; 