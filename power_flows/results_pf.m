function results_pf(pfsettings, results, num, bus, branchi, branchj, baseMVA)
fprintf('\tTOOLBOX FOR POWER SYSTEM COMPUTATIONS - POWER FLOW ANALYSIS\n')
fprintf(['\tDate: ', datestr(now, 'dd.mm.yyyy HH:MM:SS \n\n')])
fprintf('\tMethod: %s\n', results.method);
fprintf('\tPower System: %s\n', results.sys);
if results.converged
 fprintf('\tExecution time: %.2f [ms]\n', results.algtime * 1000) 
else
 fprintf('\tThe algortihm did not converged after: %d iterations!\n \t :( :( :(\n', pfsettings.maxNumberOfIter)
 return;
end
fprintf('\tNumber of iterations: %d\n', results.iter)     
if pfsettings.showbus
    disp(' ')
    fprintf('\t')   
    A = [ bus', results.Vm, results.Va .* 180/pi] ;
    disp(' ')
    disp('   _____________________________________________________________________________________')
    disp('  |                               POWER FLOW BUS DATA                                   |')
    disp('  |                                                                                     |')
    disp('  |                  VOLTAGE                 GENERATION                DEMAND           |')
    disp('  |   No.      Mag[p.u.] | Angle[deg.]    P [MW]  |  Q [MVAr]    P [MW]  |  Q [MVAr]    |')
    disp('  |  ----     ----------   -----------    -------   ---------    -------    --------    |')
    for i = 1:num.bus
        if i == num.islack
            fprintf('  |   %-6d    %6.3f      %7.2f*  ', A(i, :));
        else 
            fprintf('  |   %-6d    %6.3f      %7.2f   ', A(i, :));
        end
        if results.Pgen(i) == 0 && results.Qgen(i) == 0
            fprintf('        -         -    ');
        elseif results.Pgen(i) == 0 
            fprintf('        -     %7.2f  ', results.Qgen(i) .* baseMVA)
        elseif results.Qgen(i) == 0
            fprintf('    %7.2f       -    ', results.Pgen(i) .* baseMVA)
        else
            fprintf('    %7.2f   %7.2f  ', results.Pgen(i) * baseMVA, results.Qgen(i) .* baseMVA)
        end
        if results.Pload(i) == 0 && results.Qload(i) == 0
            fprintf('        -         -    ');
        elseif results.Pload(i) == 0 
            fprintf('        -     %7.2f  ', results.Qload(i) .* baseMVA)
        elseif results.Qload(i) == 0
            fprintf('    %7.2f       -    ', results.Pload(i) .* baseMVA)
        else
            fprintf('    %7.2f   %7.2f  ', results.Pload(i) * baseMVA, results.Qload(i) * baseMVA)
        end    
        fprintf('    |\n')                                                                
    end
    disp('  |-------------------------------------------------------------------------------------|')   
end
if pfsettings.showbranch
    fprintf('\t') 
    disp(' ')
    disp('   __________________________________________________________________________________________')
    disp('  |                               POWER FLOW BRANCH DATA                                     |')
    disp('  |                                                                                          |')
    disp('  | Branch     From   To     From Bus Power Flow    To Bus Power Flow       Power Losses     |')
    disp('  |   No.      Bus    Bus    P [MW]   | Q [MVAr]   P [MW]   | Q [MVAr]   P [MW]   | Q [MVAr] |')
    disp('  |  -----    -----  -----   --------   --------   --------   --------   --------   -------- |')
    for i = 1:num.branch
        fprintf("  |   %-6d %6d %6d", i, branchi(i), branchj(i));
        fprintf("   %8.2f   %8.2f", results.Pij(i) * baseMVA, results.Qij(i)* baseMVA);
        fprintf("   %8.2f   %8.2f", results.Pji(i) * baseMVA, results.Qji(i) * baseMVA);
        fprintf("  %7.2f   %8.2f   |\n", results.Ploss(i) * baseMVA, results.Qloss(i) * baseMVA);
    end
    disp('  |------------------------------------------------------------------------------------------|')
end
end


