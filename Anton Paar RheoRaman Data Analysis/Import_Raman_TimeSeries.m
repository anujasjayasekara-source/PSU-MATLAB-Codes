function Raman = Import_Raman_TimeSeries()

[files,path] = uigetfile( ...
    {'*.txt','Text files (*.txt)'; ...
    '*.csv','CSV files (*.csv)'; ...
    '*.*','All files'}, ...
    'Select Raman time-series files', ...
    'MultiSelect','on');

if isequal(files,0)
    Raman = [];
    return
end

if ischar(files)
    files = {files};
end

Raman = struct([]);

for i = 1:numel(files)

    filename = fullfile(path,files{i});
    fprintf('Reading %s\n',files{i})

    % Read entire file
    lines = splitlines(fileread(filename));

    % Locate spectrum section
    specStartLine = find(strcmp(strtrim(lines),'[SPECTRUM]'),1);

    specEndLine = find(strcmp(strtrim(lines),'[ANALYSIS RESULT]'),1);

    % If no [ANALYSIS RESULT] section exists, use the end of the file
    if isempty(specEndLine)
        specEndLine = numel(lines) + 1;
    end

    % Raman shift header
    StartLine = split(strtrim(lines(specStartLine+1)), ';');

    RamanShift = str2double(StartLine(3:end));
    RamanShift = RamanShift(~isnan(RamanShift));

    nRamanShifts = length(RamanShift);

    % Remaining lines are spectra
    spectrumLines = lines(specStartLine+2 : specEndLine-1);
    spectrumLines = spectrumLines(strlength(strtrim(spectrumLines)) > 0);

    nTime = length(spectrumLines);

    Intensity = zeros(nTime,nRamanShifts);
    TimeOffset = zeros(nTime,1);
    TimeStamp = strings(nTime,1);

    for j = 1:nTime

        values = split(spectrumLines(j), ';');

        % Store timestamp and time offset
        TimeStamp(j) = values(1);
        TimeOffset(j) = round(str2double(values(2)), -3)/1000;

        % Read intensity values
        numericValues = str2double(values(3:end));

        % Remove any trailing NaNs
        numericValues = numericValues(~isnan(numericValues));

        % If there is one extra numeric value, remove it
        if numel(numericValues) == nRamanShifts + 1
            numericValues = numericValues(2:end);
        end

        if numel(numericValues) ~= nRamanShifts
            error('%s spectrum %d has %d points, expected %d',...
                files{i}, j, numel(numericValues), nRamanShifts)
        end

        Intensity(j,:) = numericValues;

    end


    % Store
    Raman(i).FileName = files{i};
    Raman(i).RamanShift = RamanShift(:);
    Raman(i).TimeOffset = TimeOffset;
    Raman(i).TimeStamp = TimeStamp;

    % Matrix format:
    % column 1 = Raman shift
    % columns 2:end = spectra at each time
    Raman(i).Data = [RamanShift(:), Intensity'];

end
fprintf('All files successfully imported')
end

