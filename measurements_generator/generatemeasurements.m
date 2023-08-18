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

trackingmode = 0;
if strcmp(mgsettings.mode, 'tracking')
    trackingmode = 1;
    calcfreq = 1 / arraygcd(1 ./ ([unique(data.pmu(:, 7)); unique(data.scada(:, 5))]));
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
        hasPMU(data.pmu(i, 1)) = 1;
        measuringTimes = ceil((mgsettings.t + eps) * data.pmu(i, 7));
        mfreqPMU = mfreqPMU + measuringTimes;
        if data.pmu(i, 2) ~= - 1
            msynPMU = msynPMU + measuringTimes * (1 + numel(data.pmucurrch{i}));
        else
            msynPMU = msynPMU + measuringTimes * (1 + nAdj(i));
        end
    end
    for i = 1:nSCADA
        mSCADA = mSCADA + ceil((mgsettings.t + eps) * data.scada(i, 5));
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
            ctheta = mod(ctheta, 2 * pi);
            if ctheta > pi
                ctheta = ctheta - 2 * pi;
            end
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
    tstamps = 1;
    f = data.fn;
    loadNo = 0;
    msynPMU = 0;
    for i = 1:nPMU
        if data.pmu ~= - 1
            msynPMU = msynPMU + 1 + numel(data.pmucurrch{i});
        else
            msynPMU = msynPMU + 1 + nAdj(i);
        end 
    end
    measurements.synpmu = zeros(msynPMU, 8);
    measurements.scada = zeros(nSCADA, 6);
end
measurements.trueVoltage = complex(zeros(data.nBuses, tstamps));
c_synpmu = 1;
c_fpmu = 1;
c_scada = 1;
for i = 1:tstamps
    if trackingmode
        % dynamic settings
        dynsettings.thetaslack = thetaslack(i);
        if ~all(f == -data.fn)
            dynsettings.f = f(i);
        else
            dynsettings.f = f;
        end
        dynsettings.loadNo = loadNo;
        if loadNo == -1
            dynsettings.load = load(i);
        elseif i / calcfreq == 1/5 * mgsettings.t
            dynsettings.load = load;
        else
            dynsettings.loadNo = 0;
        end
        results = run_power_flows(pfsettings, data, dynsettings);
        
    else
        results = run_power_flows(pfsettings, data);
    end
    results_pf(data, pfsettings, results);
    % write true values
    measurements.trueVoltage(:, i) = results.Vm .* exp(1i * results.Va);

    for j = 1:nPMU
        bus = data.pmu(j, 1);
        if trackingmode
            currFreq = f(i);
            if i == 1
                currRocof = f(i) - data.fn;
            else
                currRocof = f(i) - f(i - 1);
            end
            currRocof = currRocof * calcfreq;
        end
        if ~trackingmode || mod((i-1) / calcfreq, 1/data.pmu(j, 7)) == 0
           % fpmu
           if trackingmode
            measurements.fpmu(c_fpmu, 1) = i - 1;
            measurements.fpmu(c_fpmu, 2) = bus;
            measurements.fpmu(c_fpmu, 3) = currFreq + randn * data.pmu(j, 5);
            measurements.fpmu(c_fpmu, 4) = currFreq;
            measurements.fpmu(c_fpmu, 5) = currRocof + randn * data.pmu(j, 6);
            measurements.fpmu(c_fpmu, 6) = currRocof;
            c_fpmu = c_fpmu + 1;
           end
           
           % synpmu - voltage
           measurements.synpmu(c_synpmu, 1) = i - 1;
           measurements.synpmu(c_synpmu, 2) = bus;
           measurements.synpmu(c_synpmu, 3) = 3;
           measurements.synpmu(c_synpmu, 4) = bus;
           measurements.synpmu(c_synpmu, 5) = results.Vm(bus) * (1 + randn * data.pmu(j, 3));
           measurements.synpmu(c_synpmu, 6) = results.Va(bus)  + randn * data.pmu(j, 4) * pi / 180;
           measurements.synpmu(c_synpmu, 7) = results.Vm(bus);
           measurements.synpmu(c_synpmu, 8) = results.Va(bus);
           c_synpmu = c_synpmu + 1;
           
           % synpmu - branch current(s)
           for k = 1:numel(data.pmucurrch{j})
               if data.pmucurrch{j}(k) > 0
                   line = data.pmucurrch{j}(k);
                   Iijm = results.Iijm(line);
                   Iija = results.Iija(line);
               else
                   line = -data.pmucurrch{j}(k);
                   Iijm = results.Ijim(line);
                   Iija = results.Ijia(line);
               end
             measurements.synpmu(c_synpmu, 1) = i - 1;
             measurements.synpmu(c_synpmu, 2) = bus;
             measurements.synpmu(c_synpmu, 3) = 1;
             measurements.synpmu(c_synpmu, 4) = data.pmucurrch{j}(k);
             measurements.synpmu(c_synpmu, 5) = Iijm * (1 + randn * data.pmu(j, 3) / 100);
             measurements.synpmu(c_synpmu, 6) = Iija  + randn * data.pmu(j, 4) * pi / 180;
             measurements.synpmu(c_synpmu, 7) = Iijm;
             measurements.synpmu(c_synpmu, 8) = Iija; 
             c_synpmu = c_synpmu + 1;
           end
        end
    end
    for j = 1:nSCADA
        if ~trackingmode || mod((i-1) / calcfreq, 1/data.scada(j, 5)) == 0
            measurements.scada(c_scada, 1) = i - 1;
            measurements.scada(c_scada, 2) = data.scada(j, 1);
            measurements.scada(c_scada, 3) = data.scada(j, 2);
            measurements.scada(c_scada, 4) = data.scada(j, 3);
            
            % active power flow
            switch measurements.scada(c_scada, 3)
                case 1
                    if measurements.scada(c_scada, 4) > 0
                        measurements.scada(c_scada, 5) = results.Pij(measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                        measurements.scada(c_scada, 6) = results.Pij(measurements.scada(c_scada, 4));
                    else
                        measurements.scada(c_scada, 5) = results.Pji(-measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                        measurements.scada(c_scada, 6) = results.Pji(-measurements.scada(c_scada, 4));
                    end
                case 2
                    if measurements.scada(c_scada, 4) > 0
                        measurements.scada(c_scada, 5) = results.Qij(measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                        measurements.scada(c_scada, 6) = results.Qij(measurements.scada(c_scada, 4));
                    else
                        measurements.scada(c_scada, 5) = results.Qji(-measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                        measurements.scada(c_scada, 6) = results.Qji(-measurements.scada(c_scada, 4));
                    end
                case 3
                    measurements.scada(c_scada, 5) = results.Pi(measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                    measurements.scada(c_scada, 6) = results.Pi(measurements.scada(c_scada, 4));
                case 4
                    measurements.scada(c_scada, 5) = results.Qi(measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                    measurements.scada(c_scada, 6) = results.Qi(measurements.scada(c_scada, 4));
                case 5
                    if measurements.scada(c_scada, 4) > 0
                        measurements.scada(c_scada, 5) = results.Iijm(measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                        measurements.scada(c_scada, 6) = results.Iijm(measurements.scada(c_scada, 4));
                    else
                        measurements.scada(c_scada, 5) = results.Ijim(-measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                        measurements.scada(c_scada, 6) = results.Ijim(-measurements.scada(c_scada, 4));
                    end
                case 6
                    measurements.scada(c_scada, 5) = results.Vm(measurements.scada(c_scada, 4)) + randn * data.scada(j, 4);
                    measurements.scada(c_scada, 6) = results.Vm(measurements.scada(c_scada, 4));
            end
            c_scada = c_scada + 1;
        end
    end
end
end

