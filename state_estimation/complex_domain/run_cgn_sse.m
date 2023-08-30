function [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, sesettings)
info.method = 'Gauss-Newton in Complex Variables';
% --------------------- Setting additional variables ----------------------
iter = 1;
converged = 0;
% -------------------------------------------------------------------------

% --------------------- Initialize state variables ------------------------
if sesettings.flatStart
    x = ones( 2 * data.nBuses, 1);
else
    x = [ powsys.bus.Vmi .* (1i * powsys.bus.Vai);
          powsys.bus.Vmi .* (-1i * powsys.bus.Vai);
         ];
end
% -------------------------------------------------------------------------

% ----------------------- Vector of measurement values --------------------
z = [ meas.pmu.m(meas.pmu.Ibranch) .* exp(1i .* meas.pmu.a(meas.pmu.Ibranch));
      meas.pmu.m(meas.pmu.Vnode) .* exp(1i .* meas.pmu.a(meas.pmu.Vnode));
      meas.pmu.m(meas.pmu.Ibranch) .* exp(-1i .* meas.pmu.a(meas.pmu.Ibranch));
      meas.pmu.m(meas.pmu.Vnode) .* exp(-1i .* meas.pmu.a(meas.pmu.Vnode));
      meas.pmu.m .* exp(-1i .* meas.pmu.a);
      meas.scada.m
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
W = sparse(1:2 * meas.num.pmu + meas.num.scada, 1:2 * meas.num.pmu + meas.num.scada, ...
           [ str2double(sesettings.mweights(2)) .* ones( 2 * meas.num.pmu, 1);...
           ones(meas.num.scada, 1) ]);
% -------------------------------------------------------------------------

% ------------------------- Construct C11 matrix --------------------------
lines = meas.pmu.loc(meas.pmu.Ibranch);
a = zeros(2 * meas.num.Ibranch, 1);
from = find(lines > 0);
to = find(lines < 0);
a(from + meas.num.Ibranch) = powsys.ybus.fromfrom(lines(from));
a(from + meas.num.Ibranch) = powsys.ybus.fromfrom(lines(from));
a(to + meas.num.Ibranch) = powsys.ybus.fromfrom(lines(from));
a(to + meas.num.Ibranch) = powsys.ybus.fromfrom(lines(from));
C11 = sparse( [ reshape(repelem(1:meas.num.Ibranch, 2), [], 1), ...
      meas.num.Ibranch + (1:meas.num.Vnode)], [])
% -------------------------------------------------------------------------

nNonZeroInC = 0;
for i = 1:mPMU
     if measurements.synpmu(i, 3) == 1
           nNonZeroInC = nNonZeroInC + 2;
     elseif measurements.synpmu(i, 3) == 2
           nNonZeroInC = nNonZeroInC + nnz(data.powerSystemAC.nodalMatrix...
               (measurements.synpmu(i, 2), :)); 
     elseif measurements.synpmu(i, 3) == 3
           nNonZeroInC = nNonZeroInC + 1;
     end
end

% PMU - Complex measurements jacobian
% Computed offline - linear relation
c = zeros(mPMU, 1);
elemInC = zeros(nNonZeroInC, 1);
rowInC = zeros(nNonZeroInC, 1);
colInC = zeros(nNonZeroInC, 1);
elemCounter = 1;
for i = 1:mPMU
    if measurements.synpmu(i, 3) == 1
        [ colInC(elemCounter:elemCounter + 1), elemInC(elemCounter:...
            elemCounter + 1) ] = cJ_current_flow_phasor(measurements.synpmu(i, 4), branchi, branchj, data.powerSystemAC);
        rowInC(elemCounter:elemCounter + 1) = [i, i];
        elemCounter = elemCounter + 2;
    elseif measurements.synpmu(i, 3) == 2
        [~, cols, vals] = find(data.powerSystemAC.nodalMatrix...
            (measurements.synpmu(i, 2), :));
        nNew = numel(vals);
        colInC(elemCounter:elemCounter + nNew - 1) = cols;
        elemInC(elemCounter:elemCounter + nNew - 1) = vals;
        rowInC(elemCounter:elemCounter + nNew - 1) = repmat(i, 1, nNew);
        elemCounter = elemCounter + nNew;
    elseif measurements.synpmu(i, 3) == 3
        colInC(elemCounter) = measurements.synpmu(i, 2);
        elemInC(elemCounter) = 1;
        rowInC(elemCounter) = i;
        elemCounter = elemCounter + 1;
    end
end
C = sparse(rowInC, colInC, elemInC, mPMU, data.nBuses);
zM = sparse(mPMU, data.nBuses);

% SCADA - Complex measurement Jacobian
nNonZeroInD = 0; 
for i = 1:mSCADA
    if measurements.scada(i, 3) == 1
        nNonZeroInD = nNonZeroInD + 2;
    elseif measurements.scada(i, 3) == 2
        nNonZeroInD = nNonZeroInD + 2;
    elseif measurements.scada(i, 3) == 3
        nNonZeroInD = nNonZeroInD + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 2), :));
    elseif measurements.scada(i, 3) == 4 
        nNonZeroInD = nNonZeroInD + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 2), :));
    elseif measurements.scada(i, 3) == 5
        nNonZeroInD = nNonZeroInD + 2;
    elseif measurements.scada(i, 3) == 6
        nNonZeroInD = nNonZeroInD + 1;
    end
end

% Real valued measurements
d = zeros(mSCADA, 1);
elemInD = zeros(nNonZeroInD, 1);
rowInD = zeros(nNonZeroInD, 1);
colInD = zeros(nNonZeroInD, 1);

% Iterating:
while k < sesettings.maxNumberOfIter
    for i = 1:mPMU
       if measurements.synpmu(i, 3) == 1
           c(i) = cF_current_flow_phasor(...
               measurements.synpmu(i, 4), branchi, branchj, data.powerSystemAC, x);
       elseif measurements.synpmu(i, 3) == 2
           c(i) = sum(transpose(data.powerSystemAC.nodalMatrix(...
               measurements.synpmu(i, 2), :)) .* x(1:data.nBuses));
       elseif measurements.synpmu(i, 3) == 3
           c(i) = x(measurements.synpmu(i, 2));
       end
    end
    elemCounter = 1;
    for i = 1:mSCADA
       if measurements.scada(i, 3) == 1
           [ d(i), colInD(elemCounter:elemCounter + 1), ...
               elemInD(elemCounter:elemCounter + 1) ] = ...
               c_active_power_flow(measurements.scada(i, 4), branchi, branchj, data.powerSystemAC, x);
           rowInD(elemCounter:elemCounter + 1) = [i, i];
           elemCounter = elemCounter + 2;
       elseif measurements.scada(i, 3) == 2
           [ d(i), colInD(elemCounter:elemCounter + 1), ...
               elemInD(elemCounter:elemCounter + 1) ] = ...
               c_reactive_power_flow(measurements.scada(i, 4), branchi, branchj, data.powerSystemAC, x);
           rowInD(elemCounter:elemCounter + 1) = [i, i];
           elemCounter = elemCounter + 2;
       elseif measurements.scada(i, 3) == 3
           nNew = nnz(data.powerSystemAC.nodalMatrix(measurements.scada(i, 4), :));
           [ d(i), colInD(elemCounter:elemCounter + nNew - 1), ...
               elemInD(elemCounter:elemCounter + nNew - 1) ] = ...
               c_active_power_injection(measurements.scada(i, 4), data.powerSystemAC, x);
           rowInD(elemCounter:elemCounter + nNew - 1) = repmat(i, 1, nNew);
           elemCounter = elemCounter + nNew;
       elseif measurements.scada(i, 3) == 4
           nNew = nnz(data.powerSystemAC.nodalMatrix(measurements.scada(i, 4), :));
           [ d(i), colInD(elemCounter:elemCounter + nNew - 1), ...
               elemInD(elemCounter:elemCounter + nNew - 1) ] = ...
               c_reactive_power_injection(measurements.scada(i, 4), data.powerSystemAC, x);
           rowInD(elemCounter:elemCounter + nNew - 1) = repmat(i, 1, nNew);
           elemCounter = elemCounter + nNew;
       elseif measurements.scada(i, 3) == 5
           [ d(i), colInD(elemCounter:elemCounter + 1), ...
               elemInD(elemCounter:elemCounter + 1) ] = ...
               c_branch_current_magnitude(measurements.scada(i, 4), branchi, branchj, data.powerSystemAC, x);
           rowInD(elemCounter:elemCounter + 1) = [i, i];
           elemCounter = elemCounter + 2;
       elseif measurements.scada(i, 3) == 6
           bus = measurements.scada(i, 4);
           d(i) = sqrt(x(bus) * x(bus + data.nBuses));
           colInD(elemCounter) = measurements.scada(i, 4);
           elemInD(elemCounter) = 0.5 * exp(-1i * angle(x(bus)));
           rowInD(elemCounter) = i;
           elemCounter = elemCounter + 1;
       end
    end
    
    % Combine Jacobians
    D = sparse(rowInD, colInD, elemInD, mSCADA, data.nBuses);
    H = [  C        zM
          zM      conj(C)
           D      conj(D)];
    h = [c; conj(c); d];
   
    % Solve
    r = z - h;
    transpConjH = H';
    G = transpConjH * W * H;
    g = transpConjH * W * r;
    dx = G \ g;
    x = x + dx; 
    if ~anySCADA || max(abs(dx)) < sesettings.eps
        break
    else
        k = k + 1;
    end
end
tEst = toc;
Vc = x(powsys.bus.busnew);
% Compose results:
results.nonZerosInH = nNonZeroInC + nNonZeroInD;
results.redundancy = (2 * mPMU + mSCADA)/(2 * data.nBuses);
results.t = tEst;
results.voltage = x(1:data.nBuses);
results.iter = k;
results.sys = data.case;
if k < sesettings.maxNumberOfIter
    results.converged = 1;
else
    results.converged = 0;
end
end

