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
nph = size(measurements.pmu, 1);
nf = nph;
dataPSAT.phasorMeas = measurements.pmu(:, [ 1, 2, 3, 3, 6, 7 ]);
dataPSAT.phasorMeas(:, 1) = dataPSAT.phasorMeas(:, 1) - 1;
dataPSAT.phasorMeas(dataPSAT.phasorMeas(:, 2) == 0, :) = [];
% "From" side of brances
idx = dataPSAT.phasorMeas(:, 2) == 1 & dataPSAT.phasorMeas(:, 3) > 0;
dataPSAT.phasorMeas(idx, 3) = data.branch(dataPSAT.phasorMeas(idx, 4), 1);
% "From" side of brances
idx = dataPSAT.phasorMeas(:, 2) == 1 & dataPSAT.phasorMeas(:, 3) < 0;
dataPSAT.phasorMeas(idx, 3) = data.branch(dataPSAT.phasorMeas(idx, 4), 2);
dataPSAT.freqMeas = measurements.pmu(:, [1, 2, 3, 9]);
dataPSAT.freqMeas(:, 1) = dataPSAT.freqMeas(:, 1) - 1;
dataPSAT.freqMeas(dataPSAT.freqMeas(:, 2) == 0, :) = [];
dataPSAT.freqMeas(dataPSAT.freqMeas(:, 2) == 1, :) = [];
dataPSAT.favg = measurements.f;
dataPSAT.nT = length(measurements.f);
dataPSAT.pert.time = measurements.ptime; 
dataPSAT.freq = dataPSAT.fn';
dataPSAT.tstep = 1/100;

%------------------------------Save data-----------------------------------
home = getenv('USERPROFILE');
path = strcat(home, '\psat\WorkPSAT\data\pmureport\', name, 'Bpmu');
data = dataPSAT;
save(path, 'data')
%--------------------------------------------------------------------------
