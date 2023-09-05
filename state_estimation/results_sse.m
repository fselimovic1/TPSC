function results_sse(results, powsys, meas, sesettings)
fprintf('\tTOOLBOX FOR POWER SYSTEM COMPUTATIONS - STATIC STATE ESTIMATION\n')
fprintf(['\tDate: ', datestr(now, 'dd.mm.yyyy HH:MM:SS \n\n')])
disp(' ')
fprintf('\tMethod: %s\n', results.info.method);
fprintf('\tPower System: %s\n', results.sys);
fprintf('\tMeasurement recived from SCADA (legacy): %d\n', meas.num.scada)
fprintf('\tMeasurement recived from WAMS (PMU): %d\n', meas.num.pmu)
fprintf('\tNumber of non-zeros in H: %d\n', results.info.nonZerosInH)
fprintf('\tPower system redundancy: %.2f\n', results.info.redundancy)
if results.converged
    fprintf('\tExecution time: %.2f [ms]\n', results.algtime * 1000) 
else
    fprintf('\tThe algortihm did not converged after: %d iterations!\n \t :( :( :(\n', sesettings.maxNumberOfIter)
    return;
end
fprintf('\tNumber of iterations: %d\n', results.iter) 
fprintf('\tPerformance accuracy index: %f\n', sum(abs(results.voltage - meas.Vtrue) .^ 2));
if sesettings.showresults
    disp(' ')
    fprintf('\t')
	A = [ powsys.bus.busnew, abs(results.voltage), angle(results.voltage) * 180/pi,...
          abs(meas.Vtrue), angle(meas.Vtrue) * 180/pi ...
          abs(results.voltage - meas.Vtrue) ] ;
	disp(' ')
    disp('   ____________________________________________________________________')
    disp('  |                Voltage Phasors Estimation Results                  |')
    disp('  |                                                                    |')
    disp('  |     No.         Estimated               True             Error     |')
    disp('  |--------------------------------------------------------------------|')
    fprintf('  |     %-8d    %5.3f<%-6.2f          %5.3f<%-6.2f       %.4f    |\n', A')
    disp('  |____________________________________________________________________|')
end
end
 

