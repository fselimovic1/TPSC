function [data] = run_data_generator(name, mode, freqOfOccurSCADA, ...
    freqOfOccurPMU, PMUs, occursOfLoad, loadType, fType, fSubType, fMin, tSE)
% SCADA - legacy measurements
% Measurement type (row number) - 
%               1 - Active power flow
%               2 - Reactive power flow
%               3 - Active power injection
%               4 - Rective power injection 
%               5 - Branch current magnitude
%               6 - Bus voltage magnitude   
%  1: - Frequency of occurrance 
%  2: - Standard deviation
%  3: - Reporting freqeuency
%
template.scada = [ 
                   freqOfOccurSCADA(1) 0.01 0.5
                   freqOfOccurSCADA(2) 0.01 0.5
                   freqOfOccurSCADA(3) 0.01 0.5
                   freqOfOccurSCADA(4) 0.01 0.5
                   freqOfOccurSCADA(5) 0.01 0.5
                   freqOfOccurSCADA(6) 0.01 0.5
                  ];
% MEASUREMENT DEVICES data format - SYNHROPHASOR MEASUREMENTS - WAMS
% 1: Measurement type - 
%               1 - Branch current phasor
%               2 - Injected current phasor
%               3 - Bus voltage phasor     
%
% 1: Step of appearance
% 2: Standard deviation (magnitude) [%]
% 3: Phase angle standard deviation [degrees]
% 4: Reporting frequency 10 [Hz] - coefficient
% 5: Reporting frequency 25 [Hz] - coefficient
% 6: Reporting frequency 50 [Hz] - coefficient
% !!! 4: + 5: + 6: = 1 !!!
%
template.pmu = [
                  freqOfOccurPMU(1)  0.5  0.1   PMUs(1, :)
                  freqOfOccurPMU(2)  0.5  0.1   PMUs(2, :)
                  freqOfOccurPMU(3)  0.5  0.1   PMUs(3, :)
                 ];
             
% Load Dynamics
% 1: Step of appearance
% 2: Number of possible types
template.load = [occursOfLoad 1];
  

%---------------------------About Power System-----------------------------

% about frequency
data.fn = 50;
data.fType = fType;
data.fSubType = fSubType;
data.fMin = fMin;


data.sys = name;
data.t = tSE;
data.mode = mode;
data.loadType = loadType;

%-------------------------------------------------------------------------- 

%-------------------------Generate Measurements----------------------------
[ data ] = data_acquisition(data, mpc);
%--------------------------------------------------------------------------
end

