function measurements  = data_acquisition(data, mpc)
tSE = data.t;
fn = data.fn;
fPM = 1;
if ~isempty(data.pmu)
    fPM = max(data.pmu(:, 7));
end
buses = size(data.bus, 1);
gens = mpc.gen(:, 1);
mpopt = mpoption('out.all', 0);

if strcmp(data.mode, 'static')
    samples = 1;
    data.load = [];
    data.f = fn;
else
    samples = tSE * fPM;
end
tstep = 1/fPM;

% Preallocation
Pgen = zeros(buses, 1);
Qgen = zeros(buses, 1);
phase_shift = 0;
U_all = zeros(buses, samples);
theta_all = zeros(buses, samples);

% Meausred values
measurements.pmu = zeros(samples * size(data.pmu, 1), 6);
measurements.scada = zeros(samples * size(data.scada, 1), 5);
if strcmp(data.mode, 'tracking')
    if strcmp(data.dynamics, 'const')
        data.f = repmat(data.fVal, 1, samples);
        load = loadrandomwalk(samples, data.lPerc);
    end
%     if data.fType == 3
%         data.f = get_predefined_curve(fPM, tSE, data.pstart, fn, ...
%             data.fType, data.fSubType, data.fMin);
%         load = load_for_fcurve(data.fSubType, fPM, fn, data.f);
%     elseif data.fType == 2
%         data.f = get_predefined_curve(fPM, tSE, data.pstart, fn, ...
%             data.fType, data.fSubType, data.fMin);
%         load = load_profile(fPM, tSE);
%     elseif data.fType == 1
%         data.f = repmat(fn + 1, 1, tSE * fPM);
%         load = load_profile(fPM, tSE);
%     end

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

% NOTE: when ti = k*20 + n*10 [ms] (n in N)-> angle shift is neglected for the sake
% of simplicity
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
            fI = freqProfile(i); % mean([freqProfile(i - 1), freqProfile(i)]);
            phase_shift = phase_shift + 360 * fn/fPM + 360 * fn/fPM * (fI/fn - 1);
            mpc.branch(:, 4) = mpc.branch(:, 4).*(freqProfile(i)/freqProfile(i - 1));
            mpc.branch(:, 5) = mpc.branch(:, 5).*(freqProfile(i)/freqProfile(i - 1));
            mpc.bus(:, 6) = mpc.bus(:, 6).*(freqProfile(i)/freqProfile(i - 1));
        end
    end
    % Calculations
    result = ext2int(runpf(mpc, mpopt));
  
    % Data
    U = result.bus(:, 8);
    U_all(:, i) = U;
    theta = result.bus(:, 9);
    theta = theta + repmat(phase_shift, buses, 1);
    theta(theta > 180) = theta(theta > 180) - 360;
    theta = theta.*pi/180;
    theta_all(:, i) = theta;
    Ubus= U.*cos(theta) + 1i*U.*sin(theta);
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
    for j = 1:size(data.scada, 1)
        measurements.scada((i - 1) * size(data.scada, 1) + j, 1) = i;
        measurements.scada((i - 1) * size(data.scada, 1) + j, 2) = data.scada(j, 1);
        measurements.scada((i - 1) * size(data.scada, 1) + j, 3) = data.scada(j, 2);
        if data.scada(j, 1) == 1
            nBranch = data.scada(j, 2);
            if data.scada(j, 2) > 0
                Pbranch = result.branch(nBranch, 14)./mpc.baseMVA;
            else
                nBranch = abs(nBranch);
                Pbranch = result.branch(nBranch, 16)./mpc.baseMVA;
            end
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 4) =...
                Pbranch + (randn * data.scada(j, 3));
            % Exact values
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 5) = Pbranch;   
        elseif data.scada(j, 1) == 2
            nBranch = data.scada(j, 2);
            if data.scada(j, 2) > 0
                Qbranch = result.branch(nBranch, 15)./mpc.baseMVA;
            else
                nBranch = abs(nBranch);
                Qbranch = result.branch(nBranch, 17)./mpc.baseMVA;
            end
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 4) = ...
                Qbranch + (randn * data.scada(j, 3));
            % Exact values
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 5) = Qbranch;
        elseif data.scada(j, 1) == 3
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 4) = ...
                P(data.scada(j, 2))  + randn * data.scada(j, 3);
            % Exact values
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 5) = ...
                P(data.scada(j, 2));
        elseif data.scada(j, 1) == 4
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 4) = ...
                Q(data.scada(j, 2))  + randn * data.scada(j, 3);
            % Exact values
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 5) = ...
                Q(data.scada(j, 2));  
        elseif data.scada(j, 1) == 5
            nBranch = data.scada(j, 2);
            if data.scada(j, 2) > 0
                Ibranch = (result.branch(nBranch, 14)./mpc.baseMVA - 1i *... 
                    result.branch(nBranch, 15)./mpc.baseMVA)/(conj(Ubus(result.branch(nBranch, 1))));
            else
                nBranch = abs(nBranch);
                Ibranch = (result.branch(nBranch, 16)./mpc.baseMVA - 1i * ...
                    result.branch(nBranch, 17)./mpc.baseMVA)/(conj(Ubus(result.branch(nBranch, 2))));
            end
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 4) = ...
                abs(Ibranch) + randn * data.scada(j, 3);
            % Exact values
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 5) = abs(Ibranch);
        elseif data.scada(j, 1) == 6
            % Measured values (exact + gaussian white noise)
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 4) = ...
                U(data.scada(j, 2)) + randn * data.scada(j, 3);
            % Exact values
            measurements.scada((i - 1) * size(data.pmu, 1) + j, 5) = U(data.scada(j, 2));
        end
    end
    
 % PMUs
    for j = 1:size(data.pmu, 1)
        mPeriod = 1/data.pmu(j, 7);
        if  mod(ti, mPeriod) == 0
            measurements.pmu((i - 1) * size(data.pmu, 1) + j, 1) = i;
            measurements.pmu((i - 1) * size(data.pmu, 1) + j, 2) = data.pmu(j, 1);
            measurements.pmu((i - 1) * size(data.pmu, 1) + j, 3) = data.pmu(j, 2);
            measurements.pmu((i - 1) * size(data.pmu, 1) + j, 9) = ...
                freqProfile(i) + randn * data.pmu(j, 5);
            if i == 1
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 10) = ...
                    (freqProfile(i) - fn) * fPM + randn * data.pmu(j, 6);
            else
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 10) = ...
                    (freqProfile(i) - freqProfile(i - 1)) * fPM + randn * data.pmu(j, 6);
            end
            if data.pmu(j, 1) == 1
                nBranch = data.pmu(j, 2);
                if data.pmu(j, 2) > 0
                    Ibranch = (result.branch(nBranch, 14)./mpc.baseMVA - 1i * ...
                        result.branch(nBranch, 15)./mpc.baseMVA)/(conj(Ubus(result.branch(nBranch, 1))));
                else
                    nBranch = abs(nBranch);
                    Ibranch = (result.branch(nBranch, 16)./mpc.baseMVA - ...
                        1i * result.branch(nBranch, 17)./mpc.baseMVA)/(conj(Ubus(result.branch(nBranch, 2))));
                end
                % Measured values (exact + gaussian white noise)
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 4) = ...
                    abs(Ibranch) * (1 + randn * data.pmu(j, 3)/100);
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 5) = ...
                    angle(Ibranch) + randn * data.pmu(j, 4) * pi/180;
                % Exact values
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 6) = abs(Ibranch);
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 7) = angle(Ibranch);               
            elseif data.pmu(j, 1) == 2
                Ij = I(data.pmu(j, 2));
                % Measured values (exact + gaussian white noise)
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 4) = ...
                    abs(Ij) * (1 + randn * data.pmu(j, 3)/100);
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 5) = ...
                    angle(Ij) + randn * data.pmu(j, 4) * pi/180;
                % Exact values
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 6) = abs(Ij);
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 7) = angle(Ij);                
            else
                Vj = Ubus(data.pmu(j, 2));
                % Measured values (exact + gaussian white noise)
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 4) = ...
                    abs(Vj) * (1 + randn * data.pmu(j, 3)/100);
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 5) = ...
                    angle(Vj) + randn * data.pmu(j, 4) * pi/180;
                % Exact values
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 6) = abs(Vj);
                measurements.pmu((i - 1) * size(data.pmu, 1) + j, 7) = angle(Vj);
            end
        end
    end
end
% Delete empty rows
measurements.pmu(measurements.pmu(:, 1) == 0, :) = [];
measurements.scada(measurements.scada(:, 1) == 0, :) = [];

% Save results
measurements.exactVals = U_all .* exp(1i * theta_all);
measurements.f = freqProfile;
end

