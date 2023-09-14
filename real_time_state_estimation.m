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
rtsesettings.method = 'pgne';
rtsesettings.mweights = [ "pmuscadaratio", 100 ];
rtsesettings.fc = 10;
rtsesettings.flatStart = 1;
rtsesettings.maxNumberOfIter = 50;
rtsesettings.eps = 1e-8;
rtsesettings.initialStage = 1;
rtsesettings.realtimeplot = 1;
rtsesettings.rtpbus = 15;
rtsesettings.plotpause = 0;
%--------------------------------------------------------------------------

% -------------------------Static State Estimation-------------------------
runrtse(casename, vrs, rtsesettings);
% -------------------------------------------------------------------------
