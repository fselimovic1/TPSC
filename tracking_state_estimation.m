clc
close
clearvars

%-------------------------- Tracking State Estimation ---------------------
%--------------------------------------------------------------------------

%------------------------ Generate Path Name ------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------- Power System Case ----------------------------
casename = 'case118';
vrs = 'TB';
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
tsesettings.domain = 'complex';
tsesettings.method = 'pgne';
tsesettings.mweights = [ "pmuscadaratio", 15 ];
tsesettings.fc = 50;
tsesettings.flatStart = 1;
tsesettings.maxNumberOfIter = 50;
tsesettings.eps = 1e-8;
tsesettings.initialStage = 1;
tsesettings.realtimeplot = 1;
tsesettings.rtpbus = 109;
tsesettings.plotpause = 0;
tsesettings.plotForPaper = 1;
tsesettings.measureTime = 1;
%--------------------------------------------------------------------------

% --------------------------- Run  TSE solver -----------------------------
runtse(casename, vrs, tsesettings);
% -------------------------------------------------------------------------
