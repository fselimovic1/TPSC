clc
close
clearvars

%--------------------- Transmission System Power Flows --------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case1888rte';
%--------------------------------------------------------------------------

%---------------------- Power Flows Solver - Settings ---------------------
pfsettings.domain = 'complex';
pfsettings.method = 'cgn_pf';
pfsettings.flatStart = 1;
pfsettings.maxNumberOfIter = 20;
pfsettings.eps = 1e-6;
pfsettings.postprocess = 1;
%--------------------------------------------------------------------------

%--------------------------- Power Flows ----------------------------------
run_power_flows(casename, pfsettings);
%--------------------------------------------------------------------------
