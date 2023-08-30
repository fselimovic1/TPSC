clc
close
clearvars

%---------------------------State Estimation-------------------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------- Power System Case ----------------------------
casename = 'case39';
vrs = 'A';
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
sesettings.domain = 'complex';
sesettings.method = 'cgn_sse';
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%--------------------------Static State Estimation-------------------------
runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

