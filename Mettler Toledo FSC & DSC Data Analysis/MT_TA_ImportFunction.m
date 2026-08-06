function [MT_TA, RawData] = MT_TA_ImportFunction(colTemp, colTime, tolC, rateTol)


    if nargin < 1 || isempty(colTemp), colTemp = 4; end
    if nargin < 2 || isempty(colTime), colTime = 2; end
    if nargin < 3 || isempty(tolC),    tolC    = 5; end
    if nargin < 4 || isempty(rateTol), rateTol = 0; end

    [fn, pn] = uigetfile('*.txt','Select data files','MultiSelect','on');
    if isequal(fn,0)
        MT_TA = struct('File','', 'Mass',NaN, 'Data',{{}}, 'Iso',{{}}, ...
             'Dynamic',struct('All',{{}},'Heating',{{}},'Cooling',{{}}));
        MT_TA(1)=[]; RawData = {}; return
    end
    if ischar(fn), fn = {fn}; end

    MT_TA = repmat(struct('File','', 'Mass',NaN, 'Data',{{}}, 'Iso',{{}}, ...
               'Dynamic',struct('All',{{}},'Heating',{{}},'Cooling',{{}})), 1, numel(fn));
    RawData = cell(1,numel(fn));

    for k = 1:numel(fn)
        fpath = fullfile(pn, fn{k});
        L = read_lines(fpath);
        T = table(L, strtrim(L), 'VariableNames', {'Line','Trim'});
        RawData{k} = T;

        mass = NaN;

        massIdx = find(contains(T.Line, "mg"), 1, 'last');
        if ~isempty(massIdx)
            lineStr = char(T.Line(massIdx));

            tok = regexp(lineStr, '([\d]+\.\d{4})\s*mg', 'tokens', 'once');
            if ~isempty(tok)
                mass = str2double(tok{1});
            end
        end

        nameIdx = find(~cellfun('isempty', regexpi(cellstr(T.Trim), '^\s*Curve Name:\s*$')));
        if isempty(nameIdx)
            warning('No "Curve Name:" found in %s', fn{k});
            MT_TA(k).File = fn{k};
            MT_TA(k).Data = {}; MT_TA(k).Iso = {};
            MT_TA(k).Dynamic = struct('All',{{}},'Heating',{{}},'Cooling',{{}});
            continue
        end
        blockEnd = [nameIdx(2:end)-1; height(T)];

        % Collect segments
        segs = struct('Index',{}, 'Name',{}, 'Values',{});
        for i = 1:numel(nameIdx)
            iStart = nameIdx(i); iStop = blockEnd(i);
            nameLine = nextNonEmpty_T(T, iStart+1, iStop);
            if nameLine == "", continue; end
            tok = regexp(nameLine, '\](\d+)\[', 'tokens','once');
            if isempty(tok), continue; end
            segNum = str2double(tok{1});

            rel = find(~cellfun('isempty', regexpi(cellstr(T.Trim(iStart:iStop)), '^\s*Curve Values:\s*')), 1, 'first');
            vals = [];
            if ~isempty(rel)
                curveValLine = iStart + rel - 1;
                dataStart = findFirstNumeric_T(T, curveValLine+1, iStop);
                if ~isempty(dataStart), vals = readNumericBlock_T(T, dataStart, iStop); end
            end
            segs(end+1) = struct('Index', segNum, 'Name', char(strtrim(nameLine)), 'Values', vals); 
        end
        if ~isempty(segs)
            [~,ord] = sort([segs.Index]); segs = segs(ord);
        end
        C = num2cell(segs(:));

        % tag & split (store only Values in Dynamic.*)
        dynVals = {}; heatVals = {}; coolVals = {};
        iso = {};  % keep iso as before if you still use it

        % track the most recent ISO hold time for subsequent Heating segments
        prevIsoHoldRounded = NaN;
        rateRounded = NaN;

        for j = 1:numel(C)
            seg = C{j}; V = seg.Values;
            if isempty(V) || size(V,2) < colTemp
                C{j} = seg; continue
            end

            Tm = V(:,colTemp);
            if size(V,2) >= colTime && all(isfinite(V(:,colTime))) && any(diff(V(:,colTime))~=0)
                tt = V(:,colTime);
            else
                tt = (1:size(V,1))';
            end
            good = isfinite(Tm) & isfinite(tt); Tm = Tm(good); tt = tt(good);
            if numel(Tm) < 2, C{j} = seg; continue; end

            if (max(Tm)-min(Tm)) <= tolC
                % Isothermal segment
                seg.Type    = 'iso';
                seg.IsoTemp = median(Tm,'omitnan');

                % Hold time on time axis (tt already filtered)
                if ~isempty(tt) && any(isfinite(tt))
                    holdTime = max(tt) - min(tt);
                else
                    holdTime = NaN;
                end
                % round up to next "clean" number
                if exist('roundup','file') || exist('roundup','builtin')
                    isoHoldRounded = round(holdTime, 1);
                else
                    isoHoldRounded = ceil(holdTime - 1e-12);
                end

                seg.IsoHold = holdTime;   % optional to keep in struct
                C{j} = seg;

                % Iso cell: {segmentStruct, isoTemp, isoHoldRounded}
                iso{end+1,1} = seg;              
                iso{end,  2} = seg.IsoTemp;      
                iso{end,  3} = isoHoldRounded;   

                % >>> update tracker for the NEXT Heating segment
                prevIsoHoldRounded = isoHoldRounded;


            else
                % --- dynamic case unchanged ---
                dt = diff(tt); dT = diff(Tm);
                g = isfinite(dt) & dt~=0 & isfinite(dT);
                rate = NaN;
                if any(g)
                    rate = median(dT(g)./dt(g));   % Rate calculation
                    
                end
                seg.Type = 'dynamic';  % keep Type in Data only
                C{j} = seg;

                % rounded rate: magnitude up to next int, keep sign
                if isfinite(rate)
                    rateRounded = round(abs(rate), 2);
                else
                    rateRounded = NaN;
                end

                % store only Values matrices
                dynVals{end+1,1} = V;
                if ~isnan(rate)
                    if rate > rateTol
                        heatVals{end+1,1} = V;                 
                        heatVals{end,  2} = prevIsoHoldRounded;  
                        heatVals{end,  3} = rateRounded;         
                    elseif isfinite(rate) && rate < -rateTol
                        coolVals{end+1,1} = V;                   
                        coolVals{end,  2} = prevIsoHoldRounded;  
                        coolVals{end,  3} = rateRounded;         
                    end
                end
            end

        end

        MT_TA(k).File = fn{k};
        MT_TA(k).Data = C;                       
        MT_TA(k).Iso  = iso;                     
        MT_TA(k).Dynamic.All     = dynVals;      
        MT_TA(k).Dynamic.Heating = heatVals;     
        MT_TA(k).Dynamic.Cooling = coolVals;     
        MT_TA(k).File = fn{k};
        MT_TA(k).Mass = mass;
    end
end

% --------- helpers ---------
function L = read_lines(fpath)
    try L = readlines(fpath);
    catch
        txt = fileread(fpath);
        parts = regexp(txt, '\r\n|\n|\r', 'split');
        L = string(parts(:));
    end
    L = replace(L, sprintf('\t'), ' ');
end

function s = nextNonEmpty_T(T, a, b)
    s = ""; a = max(a,1); b = min(b,height(T));
    for r = a:b
        t = strtrim(T.Line(r));
        if strlength(t) > 0, s = t; return; end
    end
end

function idx = findFirstNumeric_T(T, a, b)
    idx = []; a = max(a,1); b = min(b,height(T));
    for r = a:b
        s = strtrim(T.Line(r));
        if s == "", continue; end
        v = str2num(s); 
        if ~isempty(v) && all(isfinite(v)), idx = r; return; end
    end
end

function M = readNumericBlock_T(T, a, b)
    rows = {}; maxw = 0; a = max(a,1); b = min(b,height(T));
    for r = a:b
        s = strtrim(T.Line(r));
        if s == "", break; end
        v = str2num(s); 
        if isempty(v) || ~all(isfinite(v)), break; end
        rows{end+1} = v; 
        maxw = max(maxw, numel(v));
    end
    if isempty(rows), M = []; return; end
    M = nan(numel(rows), maxw);
    for i = 1:numel(rows)
        v = rows{i};
        M(i,1:numel(v)) = v;
    end
end

function y = roundup(x, granularity)
    if nargin < 2 || isempty(granularity), granularity = 0.1; end
    epsFix = 1e-12;
    y = ceil( (x - epsFix) ./ granularity ) .* granularity;
end

