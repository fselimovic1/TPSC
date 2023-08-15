function [ measurements ] = generatemeasurements(mgsettings, pfsettings, data)
% number of devices
nPMU = size(data.pmu, 1);
nSCADA  = size(data.scada, 1);


% Adjacency list is needed.
if ~isfield(data, 'adj')
    data.adj = adjacencylist(data);
end
nAdj = zeros(data.nBuses, 1);
for i = 1:data.nBuses
    nAdj(i) = numel(data.adj{i});
end

% In tracking mode, calculating frequency is equal to the maximum reporting     
% frequency of devices embedded in the grid.
trackingmode = 0;
if strcmp(mgsettings.mode, 'tracking')
    trackingmode = 1;
    if nPMU
        calcfreq = max(data.pmu(:, 7)); 
    else
        calcfreq = max(data.scada(:, 5));
    end
    tstamps = mgsettings.t * calcfreq;
    if round(tstamps) ~= tstamps
        tstamps = round(tstamps); 
    else
        tstamps = tstamps + 1;
    end
    % number of measurements in the interval
    msynPMU = 0;
    mfreqPMU = 0;
    mSCADA = 0;
    hasPMU = zeros(data.nBuses, 1);
    for i = 1:nPMU
        hasPMU(data.bus(i, 1)) = 1;
        mfreqPMU = mfreqPMU + fix(tstamps * data.pmu(i, 7) / calcfreq);
        if data.pmu(i, 2) ~= - 1
            msynPMU = msynPMU + fix(tstamps * data.pmu(i, 7) / calcfreq)...
                    * (1 + data.pmu(i, 2));
        else
            msynPMU = msynPMU + fix(tstamps * data.pmu(i, 7) / calcfreq)...
                    * (1 + nAdj(i));
        end
    end
    for i = 1:nSCADA
        mSCADA = mSCADA + fix(tstamps * data.scada(i, 5) / calcfreq);
    end
    measurements.synpmu = zeros(msynPMU, 8);
    measurements.fpmu = zeros(msynPMU, 6);
    measurements.scada = zeros(mSCADA, 6);
    % power system dynamics 
    fmode = mgsettings.fdynamics(1);
    lmode = mgsettings.ldynamics(1);
    if strcmp(fmode, "random") || strcmp(lmode, "random")
        walk = randomwalk(tstamps);
    end
    
    % frequency dynamics
    if strcmp(fmode, "const")
       f = -data.fn; 
    elseif strcmp(fmode, "random")
       f = data.fn * ones(tstamps, 1) + ...
                       str2double(mgsettings.fdynamics(2)) * walk;
    elseif strcmp(extractBefore(fmode, 3), "UD")
       v = str2double(extractAfter(fmode, 2));
       f = get_predefined_curve(calcfreq, mgsettings.t, mgsettings.t * ...
                       1/5, data.fn, 3, v, str2double(mgsettings.fdynamics(2)));
    end
    if ~all(f == -data.fn)
       thetaslack = zeros(tstamps, 1);
        ctheta = 0;
        for i = 1:tstamps
            ctheta = ctheta + 2 * pi * f(i) / calcfreq;
            thetaslack(i) = ctheta;
        end
    else
        thetaslack(i) = 0;
    end
    
    % load dynamics
    if strcmp(lmode, "const")
        loadNo = 0;
    elseif strcmp(lmode, "random")
        loadNo = -1;
        load = ones(tstamps, 1) - str2double(mgsettings.ldynamics(2)) * walk;
    elseif strcmp(lmode, "loadon")
        while true
            loadNo = randi(data.nBuses);
            if data.bus(loadNo, 2) == 1
                break;
            end
        end
        load = [ str2double(mgsettings.ldynamics(2)), str2double(mgsettings.ldynamics(3)) ] / data.baseMVA;
    elseif strcmp(lmode, "loadoff")
        while true
            loadNo = randi(data.nBuses);
            if data.bus(loadNo, 2) == 1
                break;
            end
        end
        load =  -[ str2double(mgsettings.ldynamics(2)), str2double(mgsettings.ldynamics(3)) ] / data.baseMVA;
    end        
else
    msynPMU = 0;
    for i = 1:nPMU
        if data.pmu ~= - 1
            msynPMU = msynPMU + 1 + data.pmu(i, 2);
        else
            msynPMU = msynPMU + 1 + nAdj(i);
        end 
    end
    measurements.synpmu = zeros(msynPMU, 8);
    measurements.fpmu = zeros(nPMU, 6);
    measurements.scada = zeros(nSCADA, 6);
end

c_synpmu = 1;
c_fpmu = 1;
c_scada = 1;
for i = 1:tstamps
    if trackingmode
        % dynamic settings
        dynsettings.thetaslack = thetaslack(i);
        if dynsettings.f ~= -data.fn
            dynsettings.f = f(i);
        else
            dynsettings.f = f;
        end
        dynsettings.loadNo = loadNo;
        if strcmp(loadNo, "all")
            dynsettings.load = load(i);
        elseif i / calcfreq == 1/5 * data.t
            dynsettings.load = load;
        else
            dynsettings.loadNo = 0;
        end
        results = run_power_flows(pfsettings, data, dynsettings);
    else
        results = run_power_flows(pfsettings, data);
    end
    for j = 1:nPMU
        bus = data.bus(j, 1);
        if ~trackingmode || mod((i-1) / calcfreq, 1/data.pmu(j, 7))
           % synpmu
           measurements.pmu(c_synpmu, 1) = i - 1;
           measurements.pmu(c_synpmu, 2) = bus;
           measurements.pmu(c_synpmu, 3) = 3;
           measurements.pmu(c_synpmu, 4) = bus;
           measurements.pmu(c_synpmu, 5) = results.Vm(bus) * (1 + randn * data.pm(j, 3));
           measurements.pmu(c_synpmu, 6) = results.Va(bus)  + randn * data.pm(j, 4);
           measurements.pmu(c_synpmu, 7) = results.Vm(bus);
           measurements.pmu(c_synpmu, 8) = results.Va(bus);
           c_synpmu = c_synpmu + 1;
           
           % 
        end
    end
    for j = 1:nSCADA
    end
end
end

