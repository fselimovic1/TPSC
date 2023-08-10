clc
close
clearvars

%--------------------- Transmission System Power Flows --------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
pfsettings.domain = 'complex';
pfsettings.method = 'cgn_pf';
pfsettings.start = 'flat';
pfsettings.maxNumberOfIter = 20;
pfsettings.eps = 1e-6;
pfsettings.postprocess = 1;
%--------------------------------------------------------------------------

%------------------------- Load Power System ------------------------------
name = 'case1354pegase';
load(strcat('SG', name));
%--------------------------------------------------------------------------

%--------------------------- Power Flows ----------------------------------
[ results ] = run_power_flows(pfsettings, data);
%--------------------------------------------------------------------------

%----------------------- Command Line Results Print -----------------------
results_pf(data, pfsettings, results)
%--------------------------------------------------------------------------
