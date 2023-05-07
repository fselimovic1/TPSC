clc
clear

%---------------------------Generate Path Name-----------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------Generate New Data-------------------------------
%--------------------Power system with measurement deiveces----------------
name = 'case300';
version = 'Dyn';
load(strcat('SG', name, 'D_', version));
mpc = ext2int(loadcase(name));
%--------------------------------------------------------------------------

%-----------------------------State Estimation mode------------------------
% 'tracking' - measurements will be taken for a specific time period
% 'static' - measurements will be taken for a single, distinct moment in time.
data.mode = 'tracking';
%--------------------------------------------------------------------------

%---------------------------Tracking SE options----------------------------
data.t = 5;
data.dynamics = 'polyfreq';
% - Constant frequency
% -            fValue
data.fParams = 50.2;
% - Curves for polynomial frequency estimation
% -              tS    A/B   fMin
data.fParams = [ 1      1    49.5 ];
data.lPerc = 0.00001;
data.occursOfLoad = 2;

%-----------------------PMU dynamic response model-------------------------
% The modified Gauss-Markov model that accounts for both the location of 
% perturbation on PMUs' response and the time elapsed.
data.dalpha = 0.5;
data.tbeta = 2;
data.k = 10;
%--------------------------------------------------------------------------

%----------------------------Generate Data---------------------------------
[ measurements ] = data_acquisition(data, mpc);
%--------------------------------------------------------------------------

%------------------------------Save Case-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\measurements\SG', name, ...
    'M_', version);
save(path, 'measurements')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------