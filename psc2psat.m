clc
clear 
close

%-----------------Power System Computations Data -> PSAT-------------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%---------------------Load Power System & Measurements---------------------
name = 'case300';
version = 'Dyn';

load(strcat('SG', name));
load(strcat('SG', name, 'M_', version))
%--------------------------------------------------------------------------

dataPSAT.nBuses = data.nBuses;
dataPSAT.nLines = data.nBranches;
dataPSAT.fn = data.fn;
dataPSAT.adj = adjacencylist(data);
dataPSAT.synpmu = measurements.synpmu;
dataPSAT.fpmu = measurements.fpmu;
dataPSAT.favg = measurements.f;
dataPSAT.nT = length(measurements.f);
dataPSAT.pert.time = 1; 
dataPSAT.freq = dataPSAT.fn';
dataPSAT.tstep = 1/50;

%------------------------------Save data-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\psat\Work\data\pmureport\', name, 'Apmu');
data = dataPSAT;
save(path, 'data')
%--------------------------------------------------------------------------
