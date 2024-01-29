clc
clear

%---------------------------Generate Path Name-----------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------Generate New Data-------------------------------
%--------------------Power system with measurement devices-----------------
casename = 'case39';
vrs = 'TC';
%--------------------------------------------------------------------------

%---------------------------- State Estimation mode -----------------------
% 'tracking' - measurements will be taken for a specific time period
% 'static' - measurements will be taken for a single, distinct moment in time.
mgsettings.mode = 'tracking';
%--------------------------------------------------------------------------

%------------------------- Tracking SE options ----------------------------
mgsettings.t = 5;
mgsettings.fdynamics = [ "const", 50.2 ];
mgsettings.ldynamics = [ "const" ];
%--------------------------------------------------------------------------

%--------------------------- Power Flow options ---------------------------
mgsettings.pfFlatStart = 0; % recommended
mgsettings.pfMaxNumOfIter = 20;
mgsettings.pfEps = 1e-6; 
%--------------------------------------------------------------------------

%----------------------- Generate Measurements ----------------------------
runmg(casename, vrs, mgsettings);
%--------------------------------------------------------------------------