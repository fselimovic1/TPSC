clc
clear

%---------------------------Generate Path Name-----------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------


%--------------------Distribute Measurement Devices------------------------
%----------------------------Power System----------------------------------
name = 'case39';
version = 'DynA';
%--------------------------------------------------------------------------

%---------------------------Legacy measurements----------------------------
% 1: Measurement type - 
%               1 - Active power flow
%               2 - Reactive power flow
%               3 - Active power injection
%               4 - Rective power injection 
%               5 - Branch current magnitude
%               6 - Bus voltage magnitude 
%                    1   2   3   4   5   6
freqOfOccurSCADA = [ 1,  1,  0,  0,  0,  1 ];
%--------------------------------------------------------------------------

%---------------------------Phasor measurements----------------------------
% freqPMUs - Frequency of occurrance
% 1: Measurement type - 
%               1 - Branch current phasor
%               2 - Injected current phasor
%               3 - Bus voltage phasor  
%                  1  2  3
freqOfOccurPMU = [ 1, 0, 1 ];

% 1 - Share of devices from P10
% 2 - Share of devices from P25
% 3 - Share of devices from P50
% 4 - Share of devices from P100
%                      1      2     3    4
percDiffRepRates = [  0.15       0.15     0.5   0.2
                      0       0     1    0
                      0.15       0.15     0.5   0.2 ];
%--------------------------------------------------------------------------

data = distribute_devices(name, freqOfOccurSCADA, freqOfOccurPMU, percDiffRepRates);
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

      