clc
close
clearvars

%-------------------- Complex Gauss Newton SE - TEST ----------------------
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case118';
vrs = 'C2';
%--------------------------------------------------------------------------

%------------------- State estimaion solver - Settings --------------------
sesettings.domain = 'complex';
sesettings.method = 'cgn_sse';
sesettings.virtual = 0;
sesettings.mweights = [ "pmuscadaratio", 5 ];
sesettings.flatStart = 0;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%------------------------ Run state estimation ----------------------------
[ resultsCNE, measurements, ~ ] = runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------
% sesettings.virtual = 1;
sesettings.method = 'cec_sse';
%------------------------ Run state estimation ----------------------------
[ resultsCEC, ~, ~ ] = runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

% ---------------------- Magnitude estimation error -----------------------
diffCNE = abs(abs(resultsCNE.voltage) - abs(measurements.trueVoltage));
diffCEC = abs(abs(resultsCEC.voltage) - abs(measurements.trueVoltage));
% -------------------------------------------------------------------------
figure
bar(diffCNE, 'linewidth', 0.5)
hold on
bar(diffCEC, 'linewidth', 0.5)
% xlabel('Bus')
% xticklabels(categories)
% legend('CGN', 'CGNW')
% ylabel('Diff')