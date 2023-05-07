clc
clear

%---------------------------Generate Path Name-----------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------


%--------------------Distribute Measurement Devices------------------------
%----------------------------Power System----------------------------------
name = 'case300';
version = 'Dyn';
%--------------------------------------------------------------------------

%---------------------------Legacy measurements----------------------------
% 1: Measurement type - 
%               1 - Active power flow
%               2 - Reactive power flow
%               3 - Active power injection
%               4 - Rective power injection 
%               5 - Branch current magnitude
%               6 - Bus voltage magnitude 
%                     1   2   3   4   5   6
scada.freqOfOccur = [ -1,  1,  1,  1,  0,  1 ];
%                      1     2     3     4     5     6 
scada.sd =          [ 0.02  0.02  0.02  0.02  0.02  0.02  ];
scada.repRate = 1;
%--------------------------------------------------------------------------

%---------------------------Phasor measurements----------------------------
% densPMUs - Density of PMUs distribution accross the grid
% When a PMU is installed on a particular bus, we assume that it is
% measuring the bus voltage, frequency and RoCoF, and a number of current 
% channels (max) is to be chosen using the variable nCurrCh (default value is -1 - all).
% 1: Measurement type - 
%               1 - Branch current phasor
%               2 - Injected current phasor
%               3 - Bus voltage phasor  
%                  1  2  3
pmu.dens = 30;
pmu.nCurrCh = -1; 
% - Standard deviations of PMU measurements
pmu.sd = [ ...
 %         mag[%]     phase[rad]   freq[Hz]   RoCoF[%]  
            0.05         0         0.002       0.04 ];
% 1 - Share of devices from P10
% 2 - Share of devices from P25
% 3 - Share of devices from P50
% 4 - Share of devices from P100
%                      1     2     3    4
pmu.percDiffRepRates = [  0.1  0.1  0.8  0  ];
%--------------------------------------------------------------------------

data = distribute_devices(name, scada, pmu);
if ~isempty(data.pmu)
    tPmu = array2table(data.pmu);
    tPmu = sortrows(tPmu, 5, 'descend');
    data.pmu = table2array(tPmu);
end
%------------------------------Save case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\power_system_with_devices\SG', name, 'D_', version);
save(path, 'data')
%--------------------------------------------------------------------------

      