function [results] = run_cgn_sse(solver, data, measurements)
tic;

anySCADA = ~isempty(measurements.scada);
data.powerSystemAC = admittance_matrix(data);

state = ones( 2 * data.nBuses, 1);
k = 1;

% Meaurements
mPMU = size(measurements.pmu, 1);
mSCADA = size(measurements.scada, 1);
z =  [ measurements.pmu(:, 4) .* exp(1i .* measurements.pmu(:, 5)) ;
       conj(measurements.pmu(:, 4) .* exp(1i .* measurements.pmu(:, 5)))
       measurements.scada(:, 4) ];
    
W = sparse(1:2 * mPMU + mSCADA, 1:2 * mPMU + mSCADA, [ 5 .* ...
    ones( 2 * mPMU, 1); ones(mSCADA, 1) ]);


nNonZeroInC = 0;
for i = 1:mPMU
     if measurements.pmu(i, 2) == 1
           nNonZeroInC = nNonZeroInC + 2;
     elseif measurements.pmu(i, 2) == 2
           nNonZeroInC = nNonZeroInC + nnz(data.powerSystemAC.nodalMatrix...
               (measurements.pmu(i, 3), :)); 
     elseif measurements.pmu(i, 2) == 3
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
    if measurements.pmu(i, 2) == 1
        [ colInC(elemCounter:elemCounter + 1), elemInC(elemCounter:...
            elemCounter + 1) ] = cJ_current_flow_phasor(measurements.pmu(i, 3), data);
        rowInC(elemCounter:elemCounter + 1) = [i, i];
        elemCounter = elemCounter + 2;
    elseif measurements.pmu(i, 2) == 2
        [~, cols, vals] = find(data.powerSystemAC.nodalMatrix...
            (measurements.pmu(i, 3), :));
        nNew = numel(vals);
        colInC(elemCounter:elemCounter + nNew - 1) = cols;
        elemInC(elemCounter:elemCounter + nNew - 1) = vals;
        rowInC(elemCounter:elemCounter + nNew - 1) = repmat(i, 1, nNew);
        elemCounter = elemCounter + nNew;
    elseif measurements.pmu(i, 2) == 3
        colInC(elemCounter) = measurements.pmu(i, 3);
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
    if measurements.scada(i, 2) == 1
        nNonZeroInD = nNonZeroInD + 2;
    elseif measurements.scada(i, 2) == 2
        nNonZeroInD = nNonZeroInD + 2;
    elseif measurements.scada(i, 2) == 3
        nNonZeroInD = nNonZeroInD + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 3), :));
    elseif measurements.scada(i, 2) == 4 
        nNonZeroInD = nNonZeroInD + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 3), :));
    elseif measurements.scada(i, 2) == 5
        nNonZeroInD = nNonZeroInD + 2;
    elseif measurements.scada(i, 2) == 6
        nNonZeroInD = nNonZeroInD + 1;
    end
end

% Real valued measurements
d = zeros(mSCADA, 1);
elemInD = zeros(nNonZeroInD, 1);
rowInD = zeros(nNonZeroInD, 1);
colInD = zeros(nNonZeroInD, 1);

% Iterating:
while k < solver.maxNumberOfIter
    for i = 1:mPMU
       if measurements.pmu(i, 2) == 1
           c(i) = cF_current_flow_phasor(...
               measurements.pmu(i, 3), data, state);
       elseif measurements.pmu(i, 2) == 2
           c(i) = sum(transpose(data.powerSystemAC.nodalMatrix(...
               measurements.pmu(i, 3), :)) .* state(1:data.nBuses));
       elseif measurements.pmu(i, 2) == 3
           c(i) = state(measurements.pmu(i, 3));
       end
    end
    elemCounter = 1;
    for i = 1:mSCADA
       if measurements.scada(i, 2) == 1
           [ d(i), colInD(elemCounter:elemCounter + 1), ...
               elemInD(elemCounter:elemCounter + 1) ] = ...
               c_active_power_flow(measurements.scada(i, 3), data, state);
           rowInD(elemCounter:elemCounter + 1) = [i, i];
           elemCounter = elemCounter + 2;
       elseif measurements.scada(i, 2) == 2
           [ d(i), colInD(elemCounter:elemCounter + 1), ...
               elemInD(elemCounter:elemCounter + 1) ] = ...
               c_reactive_power_flow(measurements.scada(i, 3), data, state);
           rowInD(elemCounter:elemCounter + 1) = [i, i];
           elemCounter = elemCounter + 2;
       elseif measurements.scada(i, 2) == 3
           nNew = nnz(data.powerSystemAC.nodalMatrix(measurements.scada(i, 3), :));
           [ d(i), colInD(elemCounter:elemCounter + nNew - 1), ...
               elemInD(elemCounter:elemCounter + nNew - 1) ] = ...
               c_active_power_injection(measurements.scada(i, 3), data, state);
           rowInD(elemCounter:elemCounter + nNew - 1) = repmat(i, 1, nNew);
           elemCounter = elemCounter + nNew;
       elseif measurements.scada(i, 2) == 4
           nNew = nnz(data.powerSystemAC.nodalMatrix(measurements.scada(i, 3), :));
           [ d(i), colInD(elemCounter:elemCounter + nNew - 1), ...
               elemInD(elemCounter:elemCounter + nNew - 1) ] = ...
               c_reactive_power_injection(measurements.scada(i, 3), data, state);
           rowInD(elemCounter:elemCounter + nNew - 1) = repmat(i, 1, nNew);
           elemCounter = elemCounter + nNew;
       elseif measurements.scada(i, 2) == 5
           [ d(i), colInD(elemCounter:elemCounter + 1), ...
               elemInD(elemCounter:elemCounter + 1) ] = ...
               c_branch_current_magnitude(measurements.scada(i, 3), data, state);
           rowInD(elemCounter:elemCounter + 1) = [i, i];
           elemCounter = elemCounter + 2;
       elseif measurements.scada(i, 2) == 6
           bus = measurements.scada(i, 3);
           d(i) = sqrt(state(bus) * state(bus + data.nBuses));
           colInD(elemCounter) = measurements.scada(i, 3);
           elemInD(elemCounter) = 0.5 * exp(-1i * angle(state(bus)));
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
    state = state + dx; 
    if ~anySCADA || max(abs(dx)) < solver.eps
        break
    else
        k = k + 1;
    end
end
tEst = toc;

% Compose results:
results.method = 'Gauss-Newton in Complex Variables';
results.nonZerosInH = nNonZeroInC + nNonZeroInD;
results.redundancy = (2 * mPMU + mSCADA)/(2 * data.nBuses);
results.t = tEst;
results.voltage = state(1:data.nBuses);
results.iter = k;
results.sys = data.case;
if k < solver.maxNumberOfIter
    results.converged = 1;
else
    results.converged = 0;
end
end

