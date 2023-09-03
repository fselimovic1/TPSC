clc
close
clearvars

% I am constantly learning about GIT!

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
rtsesettings.domain = 'real';
rtsesettings.method = 'fEKFrect';
rtsesettings.mweights = [ "pmuscadaratio", 1 ];
rtsesettings.fc = 50;
rtsesettings.flatStart = 1;
rtsesettings.maxNumberOfIter = 50;
rtsesettings.eps = 1e-6;
rtsesettings.initialStage = 1;
rtsesettings.realtimeplot = 1;
rtsesettings.rtpbus = 15;
%--------------------------------------------------------------------------

%--------------------------Static State Estimation-------------------------
runrtse(casename, vrs, rtsesettings);
%--------------------------------------------------------------------------
