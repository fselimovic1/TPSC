function results_pf(solver, results)

fprintf(['\tDate: ', datestr(now, 'dd.mm.yyyy HH:MM:SS \n')])
 disp(' ')
 fprintf('\tMethod: %s\n', results.method);
 fprintf('\tPower System: %s\n', results.sys);
 disp(' ')
 if results.converged
     fprintf('\tExecution time: %.2f [ms]\n', results.at * 1000) 
 else
     fprintf('\tThe algortihm did not converged after: %d iterations!\n \t :( :( :(\n', solver.maxNumberOfIter)
     return;
 end
 fprintf('\tNumber of iterations: %d\n', results.iter) 
 disp(' ')
 fprintf('\t')
 toShowBase = input('Show main bus data (y/n)?', 's');
 if strcmp(toShowBase, 'y') || strcmp(toShowBase, 'Y')
     nBuses = numel(results.Vm);
     A = [(1:nBuses)', results.Vm, results.Va .* 180/pi,...
           results.Pi, results.Qi, results.Pi, results.Qi ] ;
     disp(' ')
	   disp('   ______________________________________________________________________________________________________')
	   disp('  |                                      POWER FLOW BUS DATA                                             |')
	   disp('  |                                                                                                      |')
       disp('  |                         VOLTAGE                    GENERATION                  DEMAND                |')
	   disp('  |     No.          Mag[p.u.] | Angle[deg.]        P [MW] | Q [MW]            P [MW] | Q [MW]           |')
	   disp('  |------------------------------------------------------------------------------------------------------|')
	fprintf('  |     %-4d          %.3f        %4.2f             %4.2f    %4.2f              %4.2f    %4.2f           |\n', A')
	   disp('  |___________________________________________________________________________________________________________|')
 end
end

