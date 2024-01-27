clc
close
clearvars

%------------------------ Tracking State Estimation -----------------------
%--------------------------------------------------------------------------

%------------------------ Generate Path Name ------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%--------------------------- Power System Case ----------------------------
casename = 'case39';
vrs = 'TC';
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
tsesettings.domain = 'real';
tsesettings.method = 'fEKFrect';
tsesettings.mweights = [ "deviceinfo", 2 ];
tsesettings.fc = 25;
tsesettings.flatStart = 0;
tsesettings.maxNumberOfIter = 50;
tsesettings.eps = 1e-8;
tsesettings.initialStage = 1;
tsesettings.virtual = 1;
tsesettings.realtimeplot = 1;
tsesettings.rtpbus = 7;
tsesettings.plotpause = 0.05;
tsesettings.plotForPaper = 1;
tsesettings.measureTime = 1;
tsesettings.timevariantR = 0;
tsesettings.NQ = 10;
tsesettings.isQconst = 0;
%--------------------------------------------------------------------------

% --------------------------- Run  TSE solver -----------------------------
runtse(casename, vrs, tsesettings);
% -------------------------------------------------------------------------
