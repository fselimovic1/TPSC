clc
clear

%---------------------------Generate Path Name-----------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------Generate New Data-------------------------------
%--------------------Power system with measurement devices-----------------
casename = 'case118';
vrs = 'A';
%--------------------------------------------------------------------------

%---------------------------- State Estimation mode -----------------------
% 'tracking' - measurements will be taken for a specific time period
% 'static' - measurements will be taken for a single, distinct moment in time.
mgsettings.mode = 'tracking';
%--------------------------------------------------------------------------

%------------------------- Tracking SE options ----------------------------
mgsettings.t = 2;
mgsettings.fdynamics = [ "random", 0.001 ];
mgsettings.ldynamics = [ "random", 0.0005 ];
%--------------------------------------------------------------------------

%--------------------------- Power Flow options ---------------------------
mgsettings.pfFlatStart = 0; % recommended
mgsettings.pfMaxNumOfIter = 20;
mgsettings.pfEps = 1e-6; 
%--------------------------------------------------------------------------

%----------------------- Generate Measurements ----------------------------
runmg(casename, vrs, mgsettings);
%--------------------------------------------------------------------------