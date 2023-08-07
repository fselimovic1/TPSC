clc
close
clearvars

%--------------------- Transmission System Power Flows --------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
solver.domain = 'complex';
solver.method = 'cgn_pf';
solver.start = '';
solver.maxNumberOfIter = 20;
solver.eps = 1e-6;
solver.postprocess = 1;
%--------------------------------------------------------------------------

%------------------------- Load Power System ------------------------------
name = 'case9241pegase';
load(strcat('SG', name));
%--------------------------------------------------------------------------

%--------------------------- Power Flows ----------------------------------
[ results ] = run_power_flows(solver, data);
%--------------------------------------------------------------------------

%----------------------- Command Line Results Print -----------------------
results_pf(data, solver, results)
%--------------------------------------------------------------------------
