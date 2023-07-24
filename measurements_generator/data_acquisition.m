function measurements  = data_acquisition(data, mpc)
tSE = data.t;
fn = data.fn;
fPM = 1;
if ~isempty(data.pmu)
    fPM = max(data.pmu(:, 7));
end
nPMUs = size(data.pmu, 1);
nSCADA = size(data.scada, 1);
buses = size(data.bus, 1);
gens = mpc.gen(:, 1);
mpopt = mpoption('out.all', 0);

if strcmp(data.mode, 'static')
    samples = 1;
    data.load = [];
    data.f = fn;
else
    samples = tSE * fPM + 1;
end
tstep = 1/fPM;

% Preallocation
Pgen = zeros(buses, 1);
Qgen = zeros(buses, 1);
phase_shift = 0;
U_all = zeros(buses, samples);
theta_all = zeros(buses, samples);
fBus = ones(buses, 4);
pmuidx = data.pmu(:, 1);

% Meausred values
ifPMU = 1;
isynPMU = 1;
measurements.scada = zeros(samples * nSCADA, 5);

if strcmp(data.mode, 'tracking')
    if strcmp(data.dynamics, 'const')
        data.f = repmat(data.fParams, 1, samples);
        load = loadrandomwalk(samples, data.lPerc);
    elseif strcmp(data.dynamics, 'polyfreq')
        tS = data.fParams(1);
        fType = data.fParams(2);
        fMin = data.fParams(3);
        data.f = get_predefined_curve(fPM, tSE, tS, fn, 3, fType, fMin);
        load = loadrandomwalk(samples, data.lPerc);
    end
    
%   Distrubute Loads & Allocate memory
    loadBuses = find(mpc.bus(:, 3) ~= 0);
    dybuses = loadBuses(mod(1:numel(loadBuses), data.occursOfLoad) ~= 0);
    loads = numel(dybuses);
    active_dy_load = zeros(loads, samples);
    reactive_dy_load = zeros(loads, samples);
    P_load = mpc.bus(:, 3);
    if ~isempty(dybuses)
        for i=1:loads
            active_dy_load(i, :) = mpc.bus(dybuses(i), 3) * load;
            reactive_dy_load(i, :) = mpc.bus(dybuses(i), 4) * load;
        end
    end
end
% Adjacency list is needed.
if ~isfield(data, 'adj')
    data.adj = adjacencylist(data);
end

freqProfile = data.f;
for i = 1:samples
    if strcmp(data.mode, 'tracking')
        % Dynamic Load
        mpc.bus(dybuses, 3) = active_dy_load(:, i);
        mpc.bus(dybuses, 4) = reactive_dy_load(:, i);
        P_load(dybuses) = active_dy_load(:, i);
        % Frequnecy deviations
        % Model parametars
        if i == 1
            mpc.branch(:, 4) = mpc.branch(:, 4).*(freqProfile(i)/fn);
            mpc.branch(:, 5) = mpc.branch(:, 5).*(freqProfile(i)/fn);
            mpc.bus(:, 6) = mpc.bus(:, 6) .* (freqProfile(i)/fn);
        else
            fI = freqProfile(i);
            phase_shift = phase_shift + 360 * fn/fPM + 360 * fn/fPM * (fI/fn - 1);
            mpc.branch(:, 4) = mpc.branch(:, 4) .* (freqProfile(i)/freqProfile(i - 1));
            mpc.branch(:, 5) = mpc.branch(:, 5) .* (freqProfile(i)/freqProfile(i - 1));
            mpc.bus(:, 6) = mpc.bus(:, 6) .* (freqProfile(i)/freqProfile(i - 1));
        end
    end
    % Calculations
    result = ext2int(runpf(mpc, mpopt));
    % Frequency and RoCoF measurements
    if i == 1
        rocof = (freqProfile(i) - fn) * fPM;
    else
        rocof = (freqProfile(i) - freqProfile(i - 1)) * fPM;
    end
    fBus(pmuidx, 1) = ones(nPMUs, 1) .* freqProfile(i) + ...
                      randn(nPMUs, 1) .* data.pmu(:, 5);
    fBus(pmuidx, 2) = ones(nPMUs, 1) .* freqProfile(i);
    fBus(pmuidx, 3) = ones(nPMUs, 1) .* rocof + ...
                      randn(nPMUs, 1) .* data.pmu(:, 6);
    fBus(pmuidx, 4) = ones(nPMUs, 1) .* rocof;
    
    % Data
    U = result.bus(:, 8);
    U_all(:, i) = U;
    theta = result.bus(:, 9);
    theta = theta + repmat(phase_shift, buses, 1);
    theta(theta > 180) = theta(theta > 180) - 360;
    theta = theta.*pi/180;
    theta_all(:, i) = theta;
    Ubus = U.*cos(theta) + 1i*U.*sin(theta);
    Pgen(gens) = result.gen(:, 2)./mpc.baseMVA;
    Qgen(gens) = result.gen(:, 3)./mpc.baseMVA;
    Pload = result.bus(:, 3)./mpc.baseMVA;
    Qload = result.bus(:, 4)./mpc.baseMVA;
    P = Pgen - Pload;
    Q = Qgen - Qload;
    S = P + 1i.*Q;
    I = conj(S)./conj(Ubus);
    
    % Current time instant
    ti = (i - 1)/(fPM);
    % Measuring 
    % SCADA
    for j = 1:nSCADA
        measurements.scada((i - 1) * nSCADA + j, 1) = i;
        measurements.scada((i - 1) * nSCADA + j, 2) = data.scada(j, 1);
        measurements.scada((i - 1) * nSCADA + j, 3) = data.scada(j, 2);
        measurements.scada((i - 1) * nSCADA + j, 4) = data.scada(j, 3);
        if data.scada(j, 2) == 1
            nBranch = data.scada(j, 3);
            if data.scada(j, 3) > 0
                Pbranch = result.branch(nBranch, 14)./mpc.baseMVA;
            else
                nBranch = abs(nBranch);
                Pbranch = result.branch(nBranch, 16)./mpc.baseMVA;
            end
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * nSCADA + j, 5) =...
                Pbranch + (randn * data.scada(j, 4));
            % Exact values
            measurements.scada((i - 1) * nSCADA + j, 6) = Pbranch;   
        elseif data.scada(j, 2) == 2
            nBranch = data.scada(j, 3);
            if data.scada(j, 3) > 0
                Qbranch = result.branch(nBranch, 15)./mpc.baseMVA;
            else
                nBranch = abs(nBranch);
                Qbranch = result.branch(nBranch, 17)./mpc.baseMVA;
            end
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * nSCADA + j, 5) = ...
                Qbranch + (randn * data.scada(j, 4));
            % Exact values
            measurements.scada((i - 1) * nSCADA + j, 6) = Qbranch;
        elseif data.scada(j, 2) == 3
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * nSCADA + j, 5) = ...
                P(data.scada(j, 3))  + randn * data.scada(j, 4);
            % Exact values
            measurements.scada((i - 1) * nSCADA + j, 6) = ...
                P(data.scada(j, 3));
        elseif data.scada(j, 2) == 4
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * nSCADA + j, 5) = ...
                Q(data.scada(j, 3))  + randn * data.scada(j, 4);
            % Exact values
            measurements.scada((i - 1) * nSCADA + j, 6) = ...
                Q(data.scada(j, 3));  
        elseif data.scada(j, 2) == 5
            nBranch = data.scada(j, 3);
            if data.scada(j, 3) > 0
                Ibranch = (result.branch(nBranch, 14)./mpc.baseMVA - 1i *... 
                    result.branch(nBranch, 15)./mpc.baseMVA)/(conj(Ubus(result.branch(nBranch, 1))));
            else
                nBranch = abs(nBranch);
                Ibranch = (result.branch(nBranch, 16)./mpc.baseMVA - 1i * ...
                    result.branch(nBranch, 17)./mpc.baseMVA)/(conj(Ubus(result.branch(nBranch, 2))));
            end
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * nSCADA + j, 5) = ...
                abs(Ibranch) + randn * data.scada(j, 4);
            % Exact values
            measurements.scada((i - 1) * nSCADA + j, 6) = abs(Ibranch);
        elseif data.scada(j, 2) == 6
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * nSCADA + j, 5) = ...
                U(data.scada(j, 3)) + randn * data.scada(j, 4);
            % Exact values
            measurements.scada((i - 1) * nSCADA + j, 6) = U(data.scada(j, 3));
        end
    end
    
 % PMUs
    for j = 1:size(data.pmu, 1)
        mPeriod = 1/data.pmu(j, 7);
        if  mod(ti, mPeriod) == 0
            bus = data.pmu(j, 1);
            Vj = Ubus(bus);
            % - fPMU
            measurements.fpmu(ifPMU, [1, 2]) = [i, bus];
            measurements.fpmu(ifPMU, 3:6) = fBus(bus, :);
            ifPMU = ifPMU + 1;
            % - synPMU
            % - Bus voltage
            measurements.synpmu(isynPMU, 1) = i;
            measurements.synpmu(isynPMU, 2) = bus;
            measurements.synpmu(isynPMU, 3) = 3;
            measurements.synpmu(isynPMU, 4) = bus;
            measurements.synpmu(isynPMU, 5:8) = [ abs(Vj) * (1 + randn * data.pmu(j, 3)/100)...
                                  angle(Vj) + randn * data.pmu(j, 4) ...
                                  abs(Vj) angle(Vj)];
            isynPMU = isynPMU + 1;
            adjBuses = data.adj{bus};
            branches = [];
            for k = 1:numel(adjBuses)
                branches = [branches; find(data.branch(:, 1) == [bus ...
                    adjBuses(k)] & data.branch(:, 2) == [adjBuses(k) bus])];
            end
            iCh = 1;
            while data.pmu(j, 2) == -1 || iCh <= data.pmu(j, 2)
                % - Branch current
                if branches(iCh) <= data.nBranches
                    nBranch = branches(iCh);
                    Ibranch = (result.branch(nBranch, 14)./mpc.baseMVA - 1i * ...
                               result.branch(nBranch, 15)./mpc.baseMVA)/...
                               (conj(Ubus(bus)));
                else
                    nBranch = branches(iCh) - data.nBranches;
                    Ibranch = (result.branch(nBranch, 16)./mpc.baseMVA - ...
                               1i * result.branch(nBranch, 17)./mpc.baseMVA)...
                               /(conj(Ubus(bus)));
                    nBranch = -nBranch;
                end
                measurements.synpmu(isynPMU, 1) = i;
                measurements.synpmu(isynPMU, 2) = bus;
                measurements.synpmu(isynPMU, 3) = 1;
                measurements.synpmu(isynPMU, 4) = nBranch;
                measurements.synpmu(isynPMU, 5:8) = [ abs(Ibranch) * (1 + randn...
                                  * data.pmu(j, 3)/100) angle(Ibranch) + randn...
                                  * data.pmu(j, 4) abs(Ibranch) angle(Ibranch)];
                isynPMU = isynPMU + 1;              
                if iCh == numel(branches)
                    break
                end
                iCh = iCh + 1;
            end
        end
    end
end
% 0 - indexed time
measurements.scada(:, 1) = measurements.scada(:, 1) - 1;
measurements.synpmu(:, 1) = measurements.synpmu(:, 1) - 1;
measurements.fpmu(:, 1) = measurements.fpmu(:, 1) - 1;

% Save results
measurements.exactVals = U_all .* exp(1i * theta_all);
measurements.f = freqProfile;
measurements.t = tSE;
end

