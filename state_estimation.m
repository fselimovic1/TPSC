clc
close
clearvars

%-------------------------- Static State Estimation -----------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case9241pegase';
vrs = 'D';
%--------------------------------------------------------------------------

%------------------- State estimaion solver - Settings --------------------
sesettings.domain = 'complex';
sesettings.method = 'pgne';
sesettings.mweights = [ "pmuscadaratio", 5 ];
sesettings.virtual = 0;
sesettings.flatStart = 0;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%------------------------ Run state estimation ----------------------------
runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

