function [results] = run_gn_sse(solver, data, measurements)
tic;
anySCADA = ~isempty(measurements.scada);
anyPMU = ~isempty(measurements.synpmu);
data.powerSystemAC = admittance_matrix(data);
%  data = compute_branch_current_coeffs(data);

% Meaurements
mPMU = size(measurements.synpmu, 1);
mSCADA = size(measurements.scada, 1);
z =  [ measurements.synpmu(:, 5);
       measurements.synpmu(:, 6);
       measurements.scada(:, 5) ];
   
% Weigthing matrix
% 1 - Computation
% W = zeros(mPMU * 2 + mSCADA, 1);
% 2 - Simplification
W = [ 5 .* ones(2 * mPMU, 1), ones(2 * mSCADA, 1)];
% Count a number of non-zero elements in H, and compute measurement
% variances   

% Measurements Jacobian matrix
nNonZeroInH = 1;
for i = 1:mPMU
    if measurements.synpmu(i, 3) == 1
        % Branch current measurements
        nNonZeroInH = nNonZeroInH + 8;
    elseif measurements.synpmu(i, 3) == 2
    elseif measurements.synpmu(i, 3) == 3
        % Voltage phasor measurement    
        nNonZeroInH = nNonZeroInH + 2;
    end
end

for i = 1:mSCADA
    if measurements.scada(i, 3) == 1
        % Active power flow
        nNonZeroInH = nNonZeroInH + 4;
    elseif measurements.scada(i, 3) == 2
        % Reactive power flow
        nNonZeroInH = nNonZeroInH + 4;
    elseif measurements.scada(i, 3) == 3
        % Active current injection
        nNonZeroInH = nNonZeroInH + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 4), :));
    elseif measurements.scada(i, 3) == 4
        % Rective current injection
        nNonZeroInH = nNonZeroInH + nnz(...
            data.powerSystemAC.nodalMatrix(measurements.scada(i, 4), :));
    elseif measurements.scada(i, 3) == 5
         % nNonZeroInH = nNonZeroInH + 4;
    elseif measurements.scada(i, 3) == 6
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
        if measurements.synpmu(i, 3) == 1
            [ hp, colH, Hp ]  = 
        elseif measurements.synpmu(i, 3) == 2
        elseif measurements.synpmu(i, 3) == 3
            [ hp, colH, Hp ]  = 
        end
    end
end
results.t = toc;
end

