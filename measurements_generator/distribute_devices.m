function data = distribute_devices(name, ddsettings)
% load power system
load(strcat('TPSC', name));

%------------------------ SCADA MEASUREMENTS ------------------------------
% maximum number of measurements per type
maxPerType = [ 2 * data.nBranches, 2 * data.nBranches, data.nBuses, ...
               data.nBuses, 2 * data.nBranches, data.nBuses ];
           
nPerType = -1 * ones(6, 1);
sdPerType = -1 * ones(6, 1);
freqPerType = -1 * ones(6, 1);

scadanames = containers.Map({'Pij', 'Qij', 'Pi', 'Qi', 'Iij', 'Vi'}, ...
                            {1, 2, 3, 4, 5, 6});
% scada devices set
nSet = numel(ddsettings.scadaset);
i = 1;
perType = false;
if strcmp(ddsettings.scadaset(1), "perc") || strcmp(ddsettings.scadaset(1), "num")
    perType = true;
end
while i <= nSet
    toAdd = 1;
    if strcmp(ddsettings.scadaset(i), "complete")
        for j = 1:numel(nPerType)
            if nPerType(j) == -1
                nPerType(j) = maxPerType(j);
            end
        end
    end
    if ismember(ddsettings.scadaset(i), scadanames.keys)
        if ~perType
            ME = MException('MyComponent:notDefinedBehaviour', ...
                            ['It is not specified how to count SCADA devices', ...
                            ' for the specified measurement type.']);
            throw(ME);
        end
        nType = scadanames(ddsettings.scadaset(i));
        val = str2double(ddsettings.scadaset(i + 1));
        if strcmp(ddsettings.scadaset(1), "perc")
            if val < 0 || val > 100
               ME = MException('MyComponent:notDefinedBehaviour', ...
                            'A percentage must be in the range from 0 to 100.');
               throw(ME); 
            end
            nPerType(nType) = fix(val * maxPerType(nType) / 100);
        else
            if val < 0 || val > maxPerType(nType)
               ME = MException('MyComponent:notDefinedBehaviour', ...
                            'Number of devices cannot be less than zero or bigger than the maximum possible number.');
               throw(ME); 
            end
            nPerType(nType) = str2double(ddsettings.scadaset(i + 1));
        end
        toAdd = toAdd + 1;
    end
    i = i + toAdd;
end
nPerType(nPerType < 0) = 0;
% scada devices standard deviations
i = 2;
nSD = numel(ddsettings.scadasd);
isRand = 0;
isFixed = 0;
if strcmp(ddsettings.scadasd(1), "rand")
    sdPerType = -1 * ones(6, 2);
    isRand = 1;
elseif strcmp(ddsettings.scadasd(1), "fixed")
    sdPerType = -1 * ones(6, 1);
    isFixed = 1;
else
    ME = MException('MyComponent:notDefinedBehaviour', ...
                  'It is not specified a type of the standard deviation input.');
    throw(ME);
end
while i <= nSD
    toAdd = 1;
    if strcmp(ddsettings.scadasd(i), "complete")
        for j = 1:size(sdPerType, 1)
            if isRand && all(sdPerType(j, :) == [ -1, -1 ])
                sdPerType(j, 1) = str2double(ddsettings.scadasd(i + 1));
                sdPerType(j, 2) = str2double(ddsettings.scadasd(i + 2));
                toAdd = toAdd + 2;
            elseif isFixed && sdPerType(j) == -1 
                sdPerType(j) = str2double(ddsettings.scadasd(i + 1));
                toAdd = toAdd + 1;
            end
        end
    end
    if ismember(ddsettings.scadasd(i), scadanames.keys)
        if isRand
            sdPerType(scadanames(ddsettings.scadasd(i)), 1) = str2double(ddsettings.scadasd(i + 1));
            sdPerType(scadanames(ddsettings.scadasd(i)), 2) = str2double(ddsettings.scadasd(i + 2));
            toAdd = toAdd + 2;
        end
        if isFixed
            sdPerType(scadanames(ddsettings.scadasd(i))) = str2double(ddsettings.scadasd(i + 1));
            toAdd = toAdd + 1;
        end
    end
    i = i + toAdd;
end
if sum(nPerType > 0 & sdPerType == -1)
   ME = MException('MyComponent:notDefinedBehaviour', ...
                   'Standard deviation has to be defined for all measurement devices (SCADA).');
   throw(ME);
end

% scada devices reporting frequencies
i = 1;
nFreq = numel(ddsettings.scadafreq);
while i <= nFreq
    toAdd = 1;
    if strcmp(ddsettings.scadafreq(i), "complete")
        for j = 1:numel(freqPerType)
            if freqPerType(j) == -1
                freqPerType(j) = str2double(ddsettings.scadafreq(i + 1));
            end
        end
        toAdd = toAdd + 1;
    end
    if ismember(ddsettings.scadafreq(i), scadanames.keys)
        freqPerType(scadanames(ddsettings.scadafreq(i))) = str2double(ddsettings.scadafreq(i + 1));
         toAdd = toAdd + 1;
    end
    i = i + toAdd;
end
if sum(nPerType > 0 & freqPerType == -1)
   ME = MException('MyComponent:notDefinedBehaviour', ...
                   'Reporting frequncy has to be defined for all measurement devices (SCADA).');
   throw(ME);
end

% building data.scada 
data.nScada = sum(nPerType);
data.scada = zeros(data.nScada, 5);
% Pij
if isRand
    devs = rand(nPerType(1), 1) * (sdPerType(1, 2) - sdPerType(1, 1)) + sdPerType(1, 1);
else
    devs = sdPerType(1) * ones(nPerType(1), 1);
end
[ busidx, branchidx ] = randombranch(data.branch, data.nBranches, nPerType(1));
data.scada(1:nPerType(1), :) = [busidx, ones(nPerType(1), 1), branchidx,...
    devs,  freqPerType(1) * ones(nPerType(1), 1)];
nextStart = 1 + nPerType(1);
% Qij
if isRand
    devs = rand(nPerType(2), 1) * (sdPerType(2, 2) - sdPerType(2, 1)) + sdPerType(2, 1);
else
    devs = sdPerType(2) * ones(nPerType(2), 1);
end
[ busidx, branchidx ] = randombranch(data.branch, data.nBranches, nPerType(2));
data.scada(nextStart:nextStart + nPerType(2) - 1, :) = [busidx, 2 * ones(nPerType(2), 1), branchidx,...
        devs,  freqPerType(2) * ones(nPerType(2), 1)];
nextStart = nextStart + nPerType(2);
% Pi
if isRand
    devs = rand(nPerType(3), 1) * (sdPerType(3, 2) - sdPerType(3, 1)) + sdPerType(3, 1);
else
    devs = sdPerType(3) * ones(nPerType(3), 1);
end
busidx = randperm(data.nBuses);
busidx = busidx(1:nPerType(3))';
data.scada(nextStart:nextStart + nPerType(3) - 1, :) = [busidx, 3 * ones(nPerType(3), 1), busidx,...
devs,  freqPerType(3) * ones(nPerType(3), 1)];
nextStart = nextStart + nPerType(3);
% Qi
if isRand
    devs = rand(nPerType(4), 1) * (sdPerType(4, 2) - sdPerType(4, 1)) + sdPerType(4, 1);
else
    devs = sdPerType(4) * ones(nPerType(4), 1);
end
busidx = randperm(data.nBuses);
busidx = busidx(1:nPerType(4))';
data.scada(nextStart:nextStart + nPerType(4) - 1, :) = [busidx, 4 * ones(nPerType(4), 1), busidx,...
devs,  freqPerType(4) * ones(nPerType(4), 1)];
nextStart = nextStart + nPerType(4);
% Iij
if isRand
    devs = rand(nPerType(5), 1) * (sdPerType(5, 2) - sdPerType(5, 1)) + sdPerType(5, 1);
else
    devs = sdPerType(5) * ones(nPerType(5), 1);
end
[ busidx, branchidx ] = randombranch(data.branch, data.nBranches, nPerType(5));
data.scada(nextStart:nextStart + nPerType(5) - 1, :) = [busidx, 5 * ones(nPerType(5), 1), branchidx,...
    devs,  freqPerType(5) * ones(nPerType(5), 1)];
nextStart = nextStart + nPerType(5);
% Vi
if isRand
    devs = rand(nPerType(6), 1) * (sdPerType(6, 2) - sdPerType(6, 1)) + sdPerType(6, 1);
else
    devs = sdPerType(6) * ones(nPerType(6), 1);
end
busidx = randperm(data.nBuses);
busidx = busidx(1:nPerType(6))';
data.scada(nextStart:nextStart + nPerType(6) - 1, :) = [busidx, 6 * ones(nPerType(6), 1), busidx,...
devs,  freqPerType(6) * ones(nPerType(6), 1)];


% ------------------------ PMU MEASUREMENTS -------------------------------
nSet = numel(ddsettings.pmuset);
i = 1;
nCurrCh = -1;
while i <= nSet
    toAdd = 1;
    if strcmp(ddsettings.pmuset(i), "optimal")
        % To be done !
    end
    if strcmp(ddsettings.pmuset(i), "num")
        data.nPmu = str2double(ddsettings.pmuset(i + 1));
        toAdd = toAdd + 1;
    end
    if strcmp(ddsettings.pmuset(i), "perc")
        data.nPmu = fix(str2double(ddsettings.pmuset(i + 1)) * data.nBuses / 100);
        toAdd = toAdd + 1;
    end
    if strcmp(ddsettings.pmuset(i), "complete")
        data.nPmu = data.nBuses;
    end
    if strcmp(ddsettings.pmuset(i), "currCh")
        nCurrCh = str2double(ddsettings.pmuset(i + 1));
        break;
    end
    i = i + toAdd;
end

pmunames = containers.Map({'magnitude', 'angle', 'frequency', 'rocof'}, ...
                            {1, 2, 3, 4});
nSD = numel(ddsettings.pmusd);                        
isRand = 0;
isFixed = 0;
if strcmp(ddsettings.pmusd(1), "rand")
    sdPerType = -1 * ones(4, 2);
    isRand = 1;
elseif strcmp(ddsettings.pmusd(1), "fix")
    sdPerType = -1 * ones(4, 1);
    isFixed = 1;
else
    ME = MException('MyComponent:notDefinedBehaviour', ...
                  'It is not specified a type of the standard deviation input.');
    throw(ME);
end
i = 2;
while i <= nSD
    toAdd = 1;
    if ismember(ddsettings.pmusd(i), pmunames.keys)
        if isRand
            sdPerType(pmunames(ddsettings.pmusd(i)), 1) = str2double(ddsettings.pmusd(i+1));
            sdPerType(pmunames(ddsettings.pmusd(i)), 2) = str2double(ddsettings.pmusd(i+2));
            toAdd = toAdd + 2;
        elseif isFixed
            sdPerType(pmunames(ddsettings.pmusd(i))) = str2double(ddsettings.pmusd(i+1));
            toAdd = toAdd + 1;
        end
    end
    i = i + toAdd;
end
if any(sdPerType == -1)
    ME = MException('MyComponent:notDefinedBehaviour', ...
                  'Standard variance has to be insert for all quantites that a PMU measures.');
    throw(ME);
end
sdPerType(3, :) = sdPerType(3, :) ./ 1000;
pmurr = containers.Map({'P100', 'P50', 'P25', 'P10'}, ...
                            {1, 2, 3, 4});
numInSetPMU = zeros(4, 1);
i = 1;
nFreq = numel(ddsettings.pmufreq);
while i <= nFreq
    toAdd = 1;
    if strcmp(ddsettings.pmufreq(i), "complete")
        numInSetPMU(pmurr(strcat("P", ddsettings.pmufreq(i + 1)))) = data.nPmu;
        toAdd = toAdd + 1;
    end
    if ismember(ddsettings.pmufreq(i), pmurr.keys)
        numInSetPMU(pmurr(ddsettings.pmufreq(i))) = round(str2double(ddsettings.pmufreq(i + 1)) * data.nPmu / 100);
         toAdd = toAdd + 1;
    end
    i = i + toAdd;
end
if sum(numInSetPMU) ~= data.nPmu && ~strcmp(ddsettings.pmufreq(1), "complete")
   ME = MException('MyComponent:notDefinedBehaviour', ...
                   'Percentage of devices in reporting rate groups must add to 100%.');
   throw(ME);
end

% building data.pmu
pmubuses = randperm(data.nBuses);
pmubuses = pmubuses(1:data.nPmu)';
pmudevs = zeros(data.nPmu, 4);
if isRand
    for i = 1:4
        pmudevs(:, i) = rand(data.nPmu, 1) * (sdPerType(i, 2) - sdPerType(i, 1)) + sdPerType(i, 1);
    end
else
    for i = 1:4
        pmudevs(:, i) = sdPerType(i, 1) * ones(data.nPmu, 1); 
    end
end
data.pmu = [ pmubuses, nCurrCh * ones(data.nPmu, 1), pmudevs, [ 100 * ones(numInSetPMU(1), 1); 50 * ones(numInSetPMU(2), 1); ...
        25 * ones(numInSetPMU(3), 1); 10 * ones(numInSetPMU(4), 1)] ];
    
    
 % make ajdacency list
data.adj = adjacencylist(data);

% pmu current channels - priority is given to connections to 
hasPmu = zeros(data.nBuses, 1);
hasPmu(pmubuses) = 1;
data.pmucurrch = cell(data.nPmu, 1);
for i = 1:data.nPmu
    connBuses = data.adj{pmubuses(i)};
    if nCurrCh ~= -1
        toThrow = numel(connBuses) - nCurrCh;
        if toThrow > 0
            j = 1;
            while j <= numel(connBuses) && toThrow
                if hasPmu(connBuses(j))
                   connBuses(j) = [];
                   j = j - 1;
                   toThrow = toThrow - 1;
                end
                j = j + 1;
            end
            if toThrow
                connBuses(end - toThrow:end) = [];
            end
        end
    end
    for j = 1:numel(connBuses)
        nLine = [];
        nLine = find(data.branch(:, 1) == pmubuses(i) & data.branch(:, 2) == connBuses(j));
        if isempty(nLine)
            nLine = -find(data.branch(:, 2) == pmubuses(i) & data.branch(:, 1) == connBuses(j));
        end
        nLine = nLine(1);
        data.pmucurrch{i} = [ data.pmucurrch{i}, nLine];
    end
end
 
