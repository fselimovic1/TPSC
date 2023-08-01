function [results] = run_power_flows(solver, data)
if strcmp('complex', solver.domain)
    results = run_cnr_pf(solver, data);
end
end
