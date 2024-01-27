clc
close
clearvars

%---------------------------- fEKFrect TEST -------------------------------
%--------------------------------------------------------------------------
domains = [ "complex", "real", "real" ];
methods = [ "ckf_dse", "fEKFrect", "qDKF", ];
% NQ = [ 10, 20, 25 ];
magnitudeOrPhase = 1;
%--------------------------- Power System Case ----------------------------
casename = 'case39';
vrs = 'TC';
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
tsesettings.mweights = [ "deviceinfo", 5 ];
tsesettings.fc = 25;
tsesettings.flatStart = 1;
tsesettings.maxNumberOfIter = 50;
tsesettings.eps = 1e-8;
tsesettings.initialStage = 1;
tsesettings.measurementWeights = 0;
tsesettings.virtual = 0;
tsesettings.realtimeplot = 0;
tsesettings.rtpbus = 109;
tsesettings.plotpause = 0.05;
tsesettings.plotForPaper = 0;
tsesettings.measureTime = 1;
tsesettings.timevariantR = 1;
tsesettings.isQconst = 0;
tsesettings.NQ = 30;
% tsesettings.NQ = [ 5, 10, 15 ];
tsesettings.excludeStart = 5;
%--------------------------------------------------------------------------
figure
title('Performance comparison')
hold on
for i = 1:numel(domains)
    tsesettings.domain = domains(i);
    tsesettings.method = methods(i);
    % tsesettings.NQ = NQ(i);
    if i == 2
        tsesettings.virtual = 1;
    end
    % --------------------------- Run  TSE solver -------------------------
    [ results, ~, ~ ] = runtse(casename, vrs, tsesettings);
    % ---------------------------------------------------------------------
    if magnitudeOrPhase == 1
        plot(results.rmse.vm(tsesettings.excludeStart:end), '-*', 'DisplayName', [char(methods(i)), ' in ', char(domains(i)), ' domain']);
        % plot(results.rmse.vm(tsesettings.excludeStart:end), '-*', 'DisplayName', [char(methods(i)), ' in ', char(domains(i)), ' domain']);
    elseif magnitudeOrPhase == 2 
        plot(results.rmse.va(tsesettings.excludeStart:end), '-*', 'DisplayName', [char(methods(i)), ' in ', char(domains(i)), ' domain']);
    end
end
legend

