clc
close
clearvars

%-------------------------- Static State Estimation -----------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case39';
vrs = 'A';
%--------------------------------------------------------------------------

%------------------- State estimaion solver - Settings --------------------
sesettings.domain = 'complex';
sesettings.method = 'cgn_sse';
sesettings.flatStart = 1;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%------------------------ Run state estimation ----------------------------
runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

