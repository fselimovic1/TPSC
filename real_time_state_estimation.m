clc
close
clearvars

%------------------------ Real Time State Estimation ----------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------- Power System Case ----------------------------
casename = 'case39';
vrs = 'A';
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
rtsesettings.domain = 'complex';
rtsesettings.method = 'fEKFrect';
rtsesettings.fc = 10;
rtsesettings.maxNumberOfIter = 50;
rtsesettings.eps = 1e-6;
rtsesettings.showresults = 0;
%--------------------------------------------------------------------------

%--------------------------Static State Estimation-------------------------
runrtse(casename, vrs, rtsesettings);
%--------------------------------------------------------------------------