clc
close
clearvars

%--------------------- Transmission System Power Flows --------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case9';
%--------------------------------------------------------------------------

%---------------------- Power Flows Solver - Settings ---------------------
pfsettings.domain = 'complex';
pfsettings.method = 'cgn_pf';
pfsettings.flatStart = 1;
pfsettings.maxNumberOfIter = 5;
pfsettings.eps = 1e-6;
pfsettings.info = 1;
pfsettings.showbus = 0;
pfsettings.showbranch = 0;
%--------------------------------------------------------------------------

%--------------------------- Power Flows ----------------------------------
runpf(casename, pfsettings);
%--------------------------------------------------------------------------
