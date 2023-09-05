function runmg(casename, vrs, mgsettings)
measurements.mode = mgsettings.mode;
% load data
data = loadcase(casename, vrs);

% number of elements in power system
num.bus = size(data.bus, 1);
num.branch = size(data.branch, 1);
num.gen = size(data.generator, 1);
num.islack = find(data.bus(:, 2) == 3);

% settings for power flows
pfsettings.domain = 'complex';
pfsettings.method = 'cgn_pf';
pfsettings.flatStart = mgsettings.pfFlatStart;
pfsettings.maxNumberOfIter = mgsettings.pfMaxNumOfIter;
pfsettings.eps = mgsettings.pfEps;
pfsettings.info = 0;
pfsettings.showbus = 0;
pfsettings.showbranch = 0;

% Adjacency list is needed.
if ~isfield(data, 'adj')
    data.adj = adjacencylist(data);
end
nAdj = zeros(num.bus, 1);
for i = 1:num.bus
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
    hasPMU = zeros(num.bus, 1);
    for i = 1:data.nPmu
        bus = data.pmu(i, 1);
        hasPMU(bus) = 1;
        measuringTimes = ceil((mgsettings.t + eps * 10) * data.pmu(i, 7));
        mfreqPMU = mfreqPMU + measuringTimes;
        if data.pmu(i, 2) ~= - 1
            msynPMU = msynPMU + measuringTimes * (1 + numel(data.pmucurrch{i}));
        else
            msynPMU = msynPMU + measuringTimes * (1 + nAdj(bus));
        end
    end
    for i = 1:data.nScada
        mSCADA = mSCADA + ceil((mgsettings.t + eps) * data.scada(i, 5));
    end
    measurements.synpmu = zeros(msynPMU, 8);
    measurements.fpmu = zeros(mfreqPMU, 6);
    measurements.scada = zeros(mSCADA, 6);
    % power system dynamics 
    fmode = mgsettings.fdynamics(1);
    lmode = mgsettings.ldynamics(1);
    if strcmp(fmode, "random") || strcmp(lmode, "random")
        walk = randomwalk(tstamps);
    end
    
    % frequency dynamics
    if strcmp(fmode, "const")
        if numel(mgsettings.fdynamics) == 1
            f = data.fn * ones(tstamps, 1); 
        else
            f = str2double(mgsettings.fdynamics(2)) * ones(tstamps, 1);
        end
    elseif strcmp(fmode, "random")
       f = data.fn * ones(tstamps, 1) + ...
                       str2double(mgsettings.fdynamics(2)) * walk;
    elseif strcmp(extractBefore(fmode, 3), "UD")
       v = str2double(extractAfter(fmode, 2));
       f = get_predefined_curve(calcfreq, mgsettings.t, mgsettings.t * ...
                       1/5, data.fn, 3, v, str2double(mgsettings.fdynamics(2)));
    end

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
    
    
    % load dynamics
    if strcmp(lmode, "const")
        loadNo = -1;
        if numel(mgsettings.ldynamics) == 1
            load = ones(tstamps, 1);             
        else
            load = str2double(mgsettings.ldynamics(2)) .* ones(tstamps, 1);
        end
    elseif strcmp(lmode, "random")
        loadNo = -1;
        load = ones(tstamps, 1) - str2double(mgsettings.ldynamics(2)) * walk;
    elseif strcmp(lmode, "loadon")
        while true
            loadNo = randi(num.bus);
            if data.bus(loadNo, 2) == 1
                break;
            end
        end
        load = [ str2double(mgsettings.ldynamics(2)), str2double(mgsettings.ldynamics(3)) ] / data.baseMVA;
    elseif strcmp(lmode, "loadoff")
        while true
            loadNo = randi(num.bus);
            if data.bus(loadNo, 2) == 1
                break;
            end
        end
        load =  -[ str2double(mgsettings.ldynamics(2)), str2double(mgsettings.ldynamics(3)) ] / data.baseMVA;
    end        
elseif strcmp(mgsettings.mode, 'static')
    tstamps = 1;
    f = data.fn;
    loadNo = 0;
    msynPMU = 0;
    for i = 1:data.nPmu
        bus = data.pmu(i, 1);
        if data.pmu(i, 2) ~= - 1
            msynPMU = msynPMU + 1 + numel(data.pmucurrch{i});
        else
            msynPMU = msynPMU + 1 + nAdj(bus);
        end 
    end
    measurements.synpmu = zeros(msynPMU, 8);
    measurements.scada = zeros(data.nScada, 6);
else
    ME = MException('MyComponent:notDefinedBehaviour', ...
                  'Mode must be static or dynamic.');
    throw(ME);
end
pftime = zeros(tstamps, 1);
pfiter = zeros(tstamps, 1);
measurements.trueVoltage = complex(zeros(num.bus, tstamps));
c_synpmu = 1;
c_fpmu = 1;
c_scada = 1;
for i = 1:tstamps
    if trackingmode
        % dynamic settings
        dynsettings.thetaslack = thetaslack(i);
        dynsettings.f = f(i);

        dynsettings.loadNo = loadNo;
        if loadNo == -1
            dynsettings.load = load(i);
        elseif i / calcfreq > 1/5 * mgsettings.t
            dynsettings.load = load;
        else
            dynsettings.loadNo = 0;
        end
        [ results, ~ ] = runpf(casename, pfsettings, dynsettings); 
    else
        [ results, ~ ] = runpf(casename, pfsettings);
    end
    pftime(i) = results.algtime;
    pfiter(i) = results.iter;
    if ~results.converged
        ME = MException('MyComponent:notDefinedBehaviour', ...
                   'Power Flow Analysis did not converged. Measurement generation must be stopped.');
        throw(ME);
    end
    % write true values
    measurements.trueVoltage(:, i) = results.Vm .* exp(1i * results.Va);

    for j = 1:data.nPmu
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
           
            % sydata.nPmu - voltage
            measurements.synpmu(c_synpmu, 1) = i - 1;
            measurements.synpmu(c_synpmu, 2) = bus;
            measurements.synpmu(c_synpmu, 3) = 3;
            measurements.synpmu(c_synpmu, 4) = bus;
            measurements.synpmu(c_synpmu, 5) = results.Vm(bus) * (1 + randn * data.pmu(j, 3)/100);
            measurements.synpmu(c_synpmu, 6) = results.Va(bus)  + randn * data.pmu(j, 4) * pi / 180;
            measurements.synpmu(c_synpmu, 7) = results.Vm(bus);
            measurements.synpmu(c_synpmu, 8) = results.Va(bus);
            c_synpmu = c_synpmu + 1;
           
            % sydata.nPmu - branch current(s)
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
    for j = 1:data.nScada
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
if strcmp(mgsettings.mode, 'tracking')
    measurements.genFreq = calcfreq;
    measurements.tstamps = tstamps;
    measurements.f = f;
end

% Print messeage to a user
% Print a message to a user
fprintf('\tTOOLBOX FOR POWER SYSTEM COMPUTATIONS - MEASUREMENT GENERATOR\n')
fprintf(['\tDate: ', datestr(now, 'dd.mm.yyyy HH:MM:SS \n\n')])
fprintf('\tMeasurements are successfully generated for the %s in %s mode.', data.case, mgsettings.mode);
if strcmp(mgsettings.mode, 'static')
    fprintf(' Power flow calculations converged after %d iterations in %.2f seconds.\n', pfiter, pftime);
    fprintf("\t%d SCADA measurements generated.\n", data.nScada);
    fprintf("\t%d PMU (synchrophasor) measurements generated.\n", msynPMU);
else
    fprintf(' Power flow calculations converged after %.0f iterations in %.2f [ms] in average.\n', mean(pfiter), mean(pftime) * 1000);
    fprintf("\tMeasurements are generated for a time interval with a duration of %.1f seconds.\n", mgsettings.t);
    fprintf("\t%d SCADA measurements generated.\n", mSCADA);
    fprintf("\t%d PMU (synchrophasor) measurements generated.\n", msynPMU);
    fprintf("\t%d PMU (frequency and rocof) measurements generated.\n", mfreqPMU);
end
%------------------------------Save Case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\TPSC\data\measurements\TPSC', casename, ...
    'M_', vrs);
save(path, 'measurements')
%--------------------------------------------------------------------------
end

