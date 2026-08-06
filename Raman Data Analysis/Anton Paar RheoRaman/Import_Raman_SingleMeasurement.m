function Raman = Import_Raman_SingleMeasurement()
%IMPORT_RAMAN Import Raman spectra from multiple .txt files.
%
%   Raman = Import_Raman()
%
%   Prompts the user to select one or more .txt files. Data are read
%   starting immediately after the line containing [SPECTRUM]. The data
%   are assumed to be semicolon-delimited.
%
%   Output:
%       Raman(i).FileName   - Name of the file
%       Raman(i).Data       - Nx2 numeric matrix

    % Select files
    [files, path] = uigetfile('*.txt', ...
        'Select Raman spectra', ...
        'MultiSelect', 'on');

    if isequal(files,0)
        Raman = [];
        return
    end

    % Ensure files is a cell array
    if ischar(files)
        files = {files};
    end

    % Preallocate structure
    Raman = struct('FileName', cell(numel(files),1), ...
                   'Data', cell(numel(files),1));

    % Loop over files
    for i = 1:numel(files)

        filename = fullfile(path, files{i});

        % Read file as text
        txt = fileread(filename);

        % Split into lines
        lines = splitlines(txt);

        % Find [SPECTRUM]
        idx = find(strcmp(strtrim(lines), '[SPECTRUM]'), 1);

        if isempty(idx)
            warning('%s does not contain [SPECTRUM]. Skipping.', files{i});
            continue
        end

        % Keep data after [SPECTRUM]
        dataLines = lines(idx+1:end);

        % Remove blank lines
        dataLines = dataLines(strlength(strtrim(dataLines)) > 0);

        % Convert to text block
        dataText = strjoin(dataLines, newline);

        % Read semicolon-delimited numeric data
        C = textscan(dataText, '%f%f', ...
            'Delimiter', ';', ...
            'CollectOutput', true);

        % Store
        Raman(i).FileName = files{i};
        Raman(i).Data = C{1};

    end
    fprintf('All files successfully imported')
end