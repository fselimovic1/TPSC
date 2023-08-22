clc
close
clearvars

%--------------------- Transmission System Power Flows --------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%---------------------- Power Flows Solver - Settings ---------------------
pfsettings.domain = 'complex';
pfsettings.method = 'cgn_pf';
pfsettings.start = '';
pfsettings.maxNumberOfIter = 20;
pfsettings.eps = 1e-6;
pfsettings.postprocess = 1;
%--------------------------------------------------------------------------

%------------------------- Load Power System ------------------------------
name = 'case1888rte';
data = loadcase(name);
%--------------------------------------------------------------------------

%--------------------------- Power Flows ----------------------------------
[ results ] = run_power_flows(pfsettings, data);
%--------------------------------------------------------------------------

%----------------------- Command Line Results Print -----------------------
results_pf(data, pfsettings, results)
%--------------------------------------------------------------------------
