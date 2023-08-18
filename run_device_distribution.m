clc
clear

%---------------------------Generate Path Name-----------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------


%--------------------Distribute Measurement Devices------------------------
% This script enables a user to employ measurement devices into existing
% power system (variable 'name'). The variable 'version' allows creation of
% multiple different measurement sets on for a power system. All settings
% are contained as field of the structure 'ddsettings'.
%----------------------------Power System----------------------------------
name = 'case39';
version = 'A';
%--------------------------------------------------------------------------

%---------------------------Legacy measurements----------------------------
% 1: Measurement type - 
%               Pij - Active power flow
%               Qij - Reactive power flow
%               Pi - Active power injection
%               Qi - Rective power injection 
%               Iij - Branch current magnitude
%               Vi - Bus voltage magnitude

% Using the field "scadaset", a user determines SCADA measurement devices
% (RTUs) which are deployed in a grid. 
ddsettings.scadaset = ["perc", "Pij", 0, "Qij", 0, "Pi", 100, "Qi", 100, "Vi", 100 ];
% Using the field "scadavar", a user determines the variances of measurement devices
% (RTUs) which are deployed in a grid.
ddsettings.scadasd = [ "fixed", "complete", 0.005 ];

ddsettings.scadafreq = [ "complete", 1 ];
%--------------------------------------------------------------------------

%---------------------------Phasor Measurements----------------------------
% Using the field "pmuset", a user determines PMU measurement devices
% (PMUs) which are installed in a grid. 
ddsettings.pmuset = [ "num", 1, "currCh" 0 ]; %[ "perc", 95, "currCh" -1];

ddsettings.pmusd = [ "fixed", "magnitude", 0.01, 0.02, "phase", 0.2, 0.3, ...
                      "frequency", 5, 5, "rocof",  0.4, 0.4 ];

ddsettings.pmufreq = [ "complete", 10 ];
%--------------------------------------------------------------------------

data = distribute_devices(name, ddsettings);
%------------------------------Save case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\power_system_with_devices\TPSC', name, 'D_', version);
save(path, 'data')
%--------------------------------------------------------------------------

      