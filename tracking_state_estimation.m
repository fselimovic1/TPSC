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
vrs = 'TA';
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
tsesettings.domain = 'real';
tsesettings.method = 'fEKFrect';
tsesettings.mweights = [ "pmuscadaratio", 15 ];
tsesettings.fc = 10;
tsesettings.flatStart = 1;
tsesettings.maxNumberOfIter = 50;
tsesettings.eps = 1e-8;
tsesettings.initialStage = 1;
tsesettings.virtual = 1;
tsesettings.realtimeplot = 1;
tsesettings.rtpbus = 1;
tsesettings.plotpause = 0.05;
tsesettings.plotForPaper = 1;
tsesettings.measureTime = 1;
%--------------------------------------------------------------------------

% --------------------------- Run  TSE solver -----------------------------
runtse(casename, vrs, tsesettings);
% -------------------------------------------------------------------------
