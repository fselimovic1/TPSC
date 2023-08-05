 clc
 close
 clearvars

%---------------------------State Estimation-------------------------------

%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
solver.domain = 'complex';
solver.method = 'sgn_sse';
solver.fc = 10;
solver.maxNumberOfIter = 20;
solver.eps = 1e-6;
%--------------------------------------------------------------------------

%---------------------Load Power System & Measurements---------------------
name = 'case118';
version = 'A';

load(strcat('SG', name, 'D_', version));
load(strcat('SG', name, 'M_', version));
%--------------------------------------------------------------------------

%--------------------------Static State Estimation-------------------------
[ results ] = run_state_estimator(solver, data, measurements);
%--------------------------------------------------------------------------

%----------------------------State Variables-------------------------------
results_sse(solver, measurements, results)
%--------------------------------------------------------------------------
