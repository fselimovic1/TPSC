clc
close
clearvars

%-------------------------- Static State Estimation -----------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case118';
vrs = 'A';
%--------------------------------------------------------------------------

%------------------- State estimaion solver - Settings --------------------
sesettings.domain = 'complex';
sesettings.method = 'ecgn_sse';
sesettings.mweights = [ "pmuscadaratio", 5 ];
sesettings.virtual = 0;
sesettings.flatStart = 1;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%------------------------ Run state estimation ----------------------------
runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

