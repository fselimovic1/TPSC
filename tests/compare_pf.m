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
solver.start = '';
solver.maxNumberOfIter = 20;
solver.eps = 1e-8;
solver.postprocess = 1;
%--------------------------------------------------------------------------

%------------------------- Load Power System ------------------------------
name = 'case9241pegase';
load(strcat('SG', name));
%--------------------------------------------------------------------------

%------------------ Complex Newton Raphson - Power Flows ------------------
[ results ] = run_power_flows(solver, data);
%--------------------------------------------------------------------------
%--------------------------- Matpower solution ----------------------------
[ mp_results] = runpf(name);
%--------------------------------------------------------------------------

%------------------------ Comparison indices ------------------------------
fprintf("Voltage magnitude deviation. Max: %f, Average: %f\n", ...
        max(abs(results.Vm - mp_results.bus(:, 8))), mean(abs(results.Vm - mp_results.bus(:, 8))));
fprintf("Voltage angle deviation. Max: %f, Average: %f\n", ...
        max(abs(results.Va * 180 / pi - mp_results.bus(:, 9))), mean(abs(results.Va * 180 / pi - mp_results.bus(:, 9))));    
%--------------------------------------------------------------------------    