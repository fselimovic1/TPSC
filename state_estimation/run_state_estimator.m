function [results] = run_state_estimator(solver, data, measurements)
if strcmp('complex', solver.domain)
    results = run_cgn_sse(solver, data, measurements);
elseif strcmp('real', solver.domain)
    if strcmp('wls_tse', solver.method)
        results = run_wls_tse_pmu(solver, data, measurements);
    end
end
end

