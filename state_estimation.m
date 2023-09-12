clc
close
clearvars

%-------------------------- Static State Estimation -----------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case300';
vrs = 'A';
%--------------------------------------------------------------------------

%------------------- State estimaion solver - Settings --------------------
sesettings.domain = 'real';
sesettings.method = 'wls_rect_sse';
sesettings.mweights = [ "pmuscadaratio", 5 ];
sesettings.flatStart = 0;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%------------------------ Run state estimation ----------------------------
runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

