function results_sse(sesettings, measurements, results)
fprintf('\tTOOLBOX FOR POWER SYSTEM COMPUTATIONS - STATIC STATE ESTIMATION\n')
fprintf(['\tDate: ', datestr(now, 'dd.mm.yyyy HH:MM:SS \n\n')])
 disp(' ')
 fprintf('\tMethod: %s\n', results.method);
 fprintf('\tPower System: %s\n', results.sys);
 fprintf('\tMeasurement recived from SCADA (legacy): %d\n', size(measurements.scada, 1))
 fprintf('\tMeasurement recived from WAMS (PMU): %d\n', size(measurements.synpmu, 1))
 fprintf('\tNumber of non-zeros in H: %d\n', results.nonZerosInH)
 fprintf('\tPower system redundancy: %f\n', results.redundancy)
 if results.converged
     fprintf('\tExecution time: %f [ms]\n', results.t*1000) 
 else
     fprintf('\tThe algortihm did not converged after: %d iterations!\n \t :( :( :(\n', sesettings.maxNumberOfIter)
     return;
 end
 fprintf('\tNumber of iterations: %d\n', results.iter) 
 fprintf('\tPerformance accuracy index: %s\n', sum(abs(results.voltage - measurements.trueVoltage) .^ 2));
 disp(' ')
 fprintf('\t')
 toShow = input('Show nodal phasor voltages (y/n)?', 's');
 if strcmp(toShow, 'y') || strcmp(toShow, 'Y')
     nBuses = size(results.voltage);
     A = [(1:nBuses)', abs(results.voltage), angle(results.voltage) * 180/pi,...
          abs(measurements.trueVoltage), angle(measurements.trueVoltage) * 180/pi ...
          abs(results.voltage - measurements.trueVoltage) ] ;
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
 

