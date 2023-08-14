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
name = 'case300';
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
ddsettings.scadaset = [ "num", "Pi", 40, "Qi", 40, "Iij", 60 ];
% Using the field "scadavar", a user determines the variances of measurement devices
% (RTUs) which are deployed in a grid.
ddsettings.scadasd = [ "rand", "Qij", 0.01, 0.02, "Vi", 0.005, 0.1, "complete" 0.01, 0.03 ];

ddsettings.scadafreq = [ "Iij", 0.5, "complete", 2 ];
%--------------------------------------------------------------------------

%---------------------------Phasor Measurements----------------------------
% Using the field "pmuset", a user determines PMU measurement devices
% (PMUs) which are installed in a grid. 
ddsettings.pmuset = [ "perc", 10, 'currCh' 2 ];

ddsettings.pmusd = [  "rand", "magnitude", 0.01, 0.02, "angle", 0.2, 0.3, ...
                     "frequency", 5, 5, "rocof",  0.4, 0.4 ];

ddsettings.pmufreq = [ "complete", 25 ];
%--------------------------------------------------------------------------

data = distribute_devices(name, ddsettings);
%------------------------------Save case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\power_system_with_devices\TPSC', name, 'D_', version);
save(path, 'data')
%--------------------------------------------------------------------------

      