clc
close
clearvars

%-------------- Compare Power Flows Solutions to MatPower -----------------
%--------------------------------------------------------------------------

%-------------------------Generate Path Name-------------------------------
 addpath(genpath(fileparts(which(mfilename))));
%--------------------------------------------------------------------------

%----------------------------About Solver----------------------------------
solver.domain = 'complex';
solver.method = 'cgn_pf';
solver.start = 'flat';
solver.maxNumberOfIter = 20;
solver.eps = 1e-8;
%--------------------------------------------------------------------------

%------------------------- Load Power System ------------------------------
name = 'case300';
load(strcat('SG', name));
%--------------------------------------------------------------------------

%------------------ Complex Newton Raphson - Power Flows ------------------
[ results ] = run_power_flows(solver, data);
Vm = abs(results.x);
Va = angle(results.x) * 180 / pi;
%--------------------------------------------------------------------------
%--------------------------- Matpower solution ----------------------------
[ mp_results] = runpf(name);
%--------------------------------------------------------------------------

%------------------------ Comparison indices ------------------------------
fprintf("Voltage magnitude deviation. Max: %f, Average: %f\n", ...
        max(abs(Vm - mp_results.bus(:, 8))), mean(abs(Vm - mp_results.bus(:, 8))));
fprintf("Voltage angle deviation. Max: %f, Average: %f\n", ...
        max(abs(Va - mp_results.bus(:, 9))), mean(abs(Va - mp_results.bus(:, 9))));    
%--------------------------------------------------------------------------    