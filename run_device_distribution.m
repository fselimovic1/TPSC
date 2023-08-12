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
name = 'case118';
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
ddsettings.scadaset = [ "num", "Vi", 50, "complete"];

% Using the field "scadavar", a user determines the variances of measurement devices
% (RTUs) which are deployed in a grid.
ddsettings.scadasd = [ "random", "complete", 0.01, 0.04, "Pij", 0.02, 0.05 ];

ddsettings.scadafreq = [ "complete", 2, "Iij", 0.5 ];
%--------------------------------------------------------------------------

%---------------------------Phasor Measurements----------------------------
% Using the field "pmuset", a user determines PMU measurement devices
% (PMUs) which are installed in a grid. 
ddsettings.pmuset = [ "density", 85, "currCh", 3 ];

ddsettings.pmusd = [  "fixed", "magnitude", 0.01, "angle", 0.1, ...
                     "frequency", 0.005, "rocof",  0.4];

ddsettings.pmufreq = [ "complete", 10 ];
%--------------------------------------------------------------------------

data = distribute_devices(name, ddsettings);
%------------------------------Save case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\power_system_with_devices\TPSC', name, 'D_', version);
save(path, 'data')
%--------------------------------------------------------------------------

      