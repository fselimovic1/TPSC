clc
clear


%---------------------------TEST OPTIONS-----------------------------------
nRuns = 100;
arrayTime = zeros(nRuns, 1);
arrayAcc = zeros(nRuns, 1);
%--------------------------------------------------------------------------

%---------------------------Generate Path Name-----------------------------
addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------
solver.domain = 'complex';
solver.maxNumberOfIter = 20;
solver.eps = 1e-6;

%--------------------Distribute Measurement Devices------------------------
%----------------------------Power System----------------------------------
name = 'case9241pegase';
mpc = ext2int(loadcase(name));
%--------------------------------------------------------------------------

freqOfOccurSCADA = [ 0,  0,  0,  0,  0,  0 ];

freqOfOccurPMU = [ 1, 0, 1 ];

percDiffRepRates = [  0    0   1     0
                      0    0   1     0 
                      0    0   1     0 ];
                  
data = distribute_devices(name, freqOfOccurSCADA, freqOfOccurPMU, percDiffRepRates);
if ~isempty(data.pmu)
    tPmu = array2table(data.pmu);
    tPmu = sortrows(tPmu, 5, 'descend');
    data.pmu = table2array(tPmu);
end
for i = 1:nRuns
%-----------------------------State Estimation mode------------------------
    data.mode = 'static';
%--------------------------------------------------------------------------

%---------------------------Tracking SE options----------------------------
    data.occursOfLoad = 5;
    data.loadType = 'gradual';

    data.fType = 3;
    data.fSubType = 1;
    data.fMin = 48;
    data.t = 5;
%--------------------------------------------------------------------------

%----------------------------Generate Data---------------------------------
    [ measurements ] = data_acquisition(data, mpc);
%--------------------------------------------------------------------------
    [ results ] = run_state_estimator(solver, data, measurements);
    arrayTime(i) = results.t;
    arrayAcc(i) = sum(abs(results.voltage - measurements.exactVals) .^ 2);
end

fprintf('\tAverage execution time: %.3f [ms]\n', mean(arrayTime)*1000)
fprintf('\tAverage accuracy index: %.3f\n', mean(arrayAcc)*1000)