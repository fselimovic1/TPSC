clc
clear

%---------------------------Generate Path Name-----------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------Generate New Data-------------------------------
%--------------------Power system with measurement deiveces----------------
name = 'case9241pegase';
version = 'A';
load(strcat('SG', name, 'D_', version));
%--------------------------------------------------------------------------

%---------------------------- State Estimation mode -----------------------
% 'tracking' - measurements will be taken for a specific time period
% 'static' - measurements will be taken for a single, distinct moment in time.
mgsettings.mode = 'static';
%--------------------------------------------------------------------------

%------------------------- Tracking SE options ----------------------------
mgsettings.t = 5;
mgsettings.fdynamics = [ "UD1", 0.001 ];
mgsettings.ldynamics = [ "random", 0.002 ];

%--------------------------- Power Flow options ---------------------------
pfsettings.domain = 'complex';
pfsettings.method = 'cgn_pf';
pfsettings.start = '';
pfsettings.maxNumberOfIter = 20;
pfsettings.eps = 1e-6;
pfsettings.postprocess = 1;
%--------------------------------------------------------------------------

%----------------------------Generate Data---------------------------------
[ measurements ] = generatemeasurements(mgsettings, pfsettings, data);
%--------------------------------------------------------------------------

%------------------------------Save Case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\measurements\SG', name, ...
    'M_', version);
save(path, 'measurements')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------