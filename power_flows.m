clc
close
clearvars

%--------------------- Transmission System Power Flows --------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case9241pegase';
%--------------------------------------------------------------------------

%---------------------- Power Flows Solver - Settings ---------------------
pfsettings.domain = 'complex';
pfsettings.method = 'cgn_pf';
pfsettings.flatStart = 0;
pfsettings.maxNumberOfIter = 20;
pfsettings.eps = 1e-6;
pfsettings.postprocess = 1;
pfsettings.showbus = 0;
pfsettings.showbranch = 0;
%--------------------------------------------------------------------------

%--------------------------- Power Flows ----------------------------------
runpf(casename, pfsettings);
%--------------------------------------------------------------------------
