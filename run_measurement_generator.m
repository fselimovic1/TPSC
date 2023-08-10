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

%-----------------------------State Estimation mode------------------------
% 'tracking' - measurements will be taken for a specific time period
% 'static' - measurements will be taken for a single, distinct moment in time.
data.mode = 'static';
%--------------------------------------------------------------------------

%------------------------- Tracking SE options ----------------------------
data.t = 5;
data.fdynamics = [ "UD1", 0.001 ];
data.ldynamics = [ "random", 0.002 ];
%----------------------------Generate Data---------------------------------
[ measurements ] = generatemeasurements(data);
%--------------------------------------------------------------------------

%------------------------------Save Case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\measurements\SG', name, ...
    'M_', version);
save(path, 'measurements')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------