clc
close
clearvars

%-------------------- Complex Gauss Newton SE - TEST ----------------------
%--------------------------------------------------------------------------

%------------------------- Power System Case ------------------------------
casename = 'case118';
vrs = 'C';
%--------------------------------------------------------------------------

%------------------- State estimaion solver - Settings --------------------
sesettings.domain = 'complex';
sesettings.method = 'cgn_sse';
sesettings.virtual = 0;
sesettings.mweights = [ "pmuscadaratio", 1 ];
sesettings.flatStart = 0;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%------------------------ Run state estimation ----------------------------
[results, measurements, ~ ] = runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

sesettings.mweights = [ "pmuscadaratio", 5 ];
%------------------------ Run state estimation ----------------------------
[resultsW, ~, ~ ] = runsse(casename, vrs, sesettings);
%--------------------------------------------------------------------------

% ---------------------- Magnitude estimation error -----------------------
diff = abs(abs(results.voltage) - abs(measurements.trueVoltage));
diffW = abs(abs(resultsW.voltage) - abs(measurements.trueVoltage));
% -------------------------------------------------------------------------
figure
bar(diffW)
hold on
bar(diff)
% xlabel('Bus')
% xticklabels(categories)
% legend('CGN', 'CGNW')
% ylabel('Diff')



