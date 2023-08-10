function results_pf(data, pfsettings, results)
fprintf('\tTOOLBOX FOR POWER SYSTEM COMPUTATIONS\n')
fprintf(['\tDate: ', datestr(now, 'dd.mm.yyyy HH:MM:SS \n\n')])
 fprintf('\tMethod: %s\n', results.method);
 fprintf('\tPower System: %s\n', results.sys);
 if results.converged
     fprintf('\tExecution time: %.2f [ms]\n', results.at * 1000) 
 else
     fprintf('\tThe algortihm did not converged after: %d iterations!\n \t :( :( :(\n', pfsettings.maxNumberOfIter)
     return;
 end
 fprintf('\tNumber of iterations: %d\n', results.iter) 
 disp(' ')
 fprintf('\t')
 toShowBase = input('Show main bus data (y/n)?', 's');
 if strcmp(toShowBase, 'y') || strcmp(toShowBase, 'Y')
     nBuses = numel(results.Vm);
     A = [ (1:nBuses)', results.Vm, results.Va .* 180/pi] ;
     disp(' ')
	   disp('   _____________________________________________________________________________________')
	   disp('  |                               POWER FLOW BUS DATA                                   |')
	   disp('  |                                                                                     |')
       disp('  |                  VOLTAGE                 GENERATION                DEMAND           |')
	   disp('  |   No.      Mag[p.u.] | Angle[deg.]    P [MW]  |  Q [MVAr]    P [MW]  |  Q [MVAr]    |')
	   disp('  |  ----     ----------   -----------    -------   ---------    -------    --------    |')
       for i = 1:nBuses
           if i == data.slackNo
               fprintf('  |   %-6d    %6.3f      %7.2f*  ', A(i, :));
           else 
               fprintf('  |   %-6d    %6.3f      %7.2f   ', A(i, :));
           end
            if results.Pgen(i) == 0 && results.Qgen(i) == 0
                fprintf('        -         -    ');
            elseif results.Pgen(i) == 0 
                fprintf('        -     %7.2f  ', results.Qgen(i) .* data.baseMVA)
            elseif results.Qgen(i) == 0
                fprintf('    %7.2f       -    ', results.Pgen(i) .* data.baseMVA)
            else
                fprintf('    %7.2f   %7.2f  ', results.Pgen(i) * data.baseMVA, results.Qgen(i) .* data.baseMVA)
            end
            if results.Pload(i) == 0 && results.Qload(i) == 0
                fprintf('        -         -    ');
            elseif results.Pload(i) == 0 
                fprintf('        -     %7.2f  ', results.Qload(i) .* data.baseMVA)
            elseif results.Qload(i) == 0
                fprintf('    %7.2f       -    ', results.Pload(i) .* data.baseMVA)
            else
                fprintf('    %7.2f   %7.2f  ', results.Pload(i) * data.baseMVA, results.Qload(i) * data.baseMVA)
            end    
            fprintf('    |\n')                                                                
        end
        disp('  |-------------------------------------------------------------------------------------|')   
 end
 fprintf('\t')
toShowFlows= input('Show power flows accross the branches (y/n)?', 's');
if strcmp(toShowFlows, 'y') || strcmp(toShowFlows, 'Y')
%     A = [ (1:data.nBranches)', data.branch(:, 1), data.powerSystemAC.to, ...
%           results.Pij,  results.Qij, results.Pji,  results.Qji, ...
%           results.Ploss,  results.Qloss ] ;
       disp(' ')
	   disp('   __________________________________________________________________________________________')
	   disp('  |                               POWER FLOW BRANCH DATA                                     |')
	   disp('  |                                                                                          |')
       disp('  | Branch     From   To     From Bus Power Flow    To Bus Power Flow       Power Losses     |')
	   disp('  |   No.      Bus    Bus    P [MW]   | Q [MVAr]   P [MW]   | Q [MVAr]   P [MW]   | Q [MVAr] |')
	   disp('  |  -----    -----  -----   --------   --------   --------   --------   --------   -------- |')
       for i = 1:data.nBranches
           fprintf("  |   %-6d %6d %6d", i, data.branch(i, 1), data.branch(i, 2));
           fprintf("   %8.2f   %8.2f", results.Pij(i) * data.baseMVA, results.Qij(i)* data.baseMVA);
           fprintf("   %8.2f   %8.2f", results.Pji(i) * data.baseMVA, results.Qji(i) * data.baseMVA);
           fprintf("  %7.2f   %8.2f   |\n", results.Ploss(i) * data.baseMVA, results.Qloss(i) * data.baseMVA);
       end
       disp('  |------------------------------------------------------------------------------------------|')
end
end

