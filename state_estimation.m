 clc
 close
 clearvars

%---------------------------State Estimation-------------------------------

%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
sesettings.domain = 'complex';
sesettings.method = 'sgn_sse';
sesettings.fc = 10;
sesettings.maxNumberOfIter = 20;
sesettings.eps = 1e-6;
%--------------------------------------------------------------------------

%---------------------Load Power System & Measurements---------------------
name = 'case118';
version = 'A';

load(strcat('SG', name, 'D_', version));
load(strcat('SG', name, 'M_', version));
%--------------------------------------------------------------------------

%--------------------------Static State Estimation-------------------------
[ results ] = run_state_estimator(sesettings, data, measurements);
%--------------------------------------------------------------------------

%----------------------------State Variables-------------------------------
results_sse(sesettings, measurements, results)
%--------------------------------------------------------------------------
