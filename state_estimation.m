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
sesettings.method = 'cgn_sse';
sesettings.fc = 10;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
%--------------------------------------------------------------------------

%---------------------Load Power System & Measurements---------------------
name = 'case39';
version = 'A';

load(strcat('TPSC', name, 'D_', version));
load(strcat('TPSC', name, 'M_', version));
%--------------------------------------------------------------------------

%--------------------------Static State Estimation-------------------------
[ results ] = run_state_estimator(sesettings, data, measurements);
%--------------------------------------------------------------------------

%----------------------------State Variables-------------------------------
results_sse(sesettings, measurements, results)
%--------------------------------------------------------------------------
