function [results] = run_gn_sse(solver, data, measurements)
tic;
anySCADA = ~isempty(measurements.scada);
anyPMU = ~isempty(measurements.pmu);
data.powerSystemAC = admittance_matrix(data);
data = compute_branch_current_coeffs(data);

% Meaurements
mPMU = size(measurements.pmu, 1);
mSCADA = size(measurements.scada, 1);
z =  [ measurements.pmu(:, 4);
       measurements.pmu(:, 5);
       measurements.scada(:, 4) ];
% Weigthing matrix

% Measurements Jacobian matrix
nNonZeroInH = 1;

for i = 1:mPMU
    if measurements.pmu(i, 2) == 1
        % Branch current measurements
        nNonZeroInH = nNonZeroInH + 8;
    elseif measurements.pmu(i, 2) == 2
    elseif measurements.pmu(i, 2) == 3
        % Voltage phasor measurement    
        nNonZeroInH = nNonZeroInH + 2;
    end
end

for i = 1:mSCADA
    if measurements.scada(i, 2) == 1
        % Active power flow
        nNonZeroInH = nNonZeroInH + 4;
    elseif measurements.scada(i, 2) == 2
        % Reactive power flow
        nNonZeroInH = nNonZeroInH + 4;
    elseif measurements.scada(i, 2) == 3
        % Active current injection
        nNonZeroInH = nNonZeroInH + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 3), :));
    elseif measurements.scada(i, 2) == 4
        % Rective current injection
        nNonZeroInH = nNonZeroInH + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 3), :));
    elseif measurements.scada(i, 2) == 5
         nNonZeroInH = nNonZeroInH + 4;
    elseif measurements.scada(i, 2) == 6
        nNonZeroInH = nNonZeroInH + 1;
    end
end
state = [ zeros(data.nBuses, 1); ones(data.nBuses, 1) ];
k = 1;
h = zeros( 2 * mPMU + mSCADA, 1);
elemInH = zeros(nNonZeroInH, 1);
rowInH = zeros(nNonZeroInH, 1);
colInH = zeros(nNonZeroInH, 1);

while k < solver.maxNumberOfIter
    for i = 1:mPMU
        if measurements.pmu(i, 2) == 1
            [ hp, colH, Hp ]  = 
        elseif measurements.pmu(i, 2) == 2
        elseif measurements.pmu(i, 2) == 3
            [ hp, colH, Hp ]  = 
        end
    end
end
results.t = toc;
end

