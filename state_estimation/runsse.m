function [ results, measurements, data ] = runsse(casename, vrs, sesettings)
% ----------------- Load Power System Data and Measurements ---------------
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');
% -------------------------------------------------------------------------

%----------------- Extract Useful Informations (Power System) -------------
powsys = preprocess_ps(data, 'se');
%--------------------------------------------------------------------------

%----------------- Extract Useful Informations (Measurements) -------------
meas = preprocess_meas(measurements);
%--------------------------------------------------------------------------

% ---------------- Solve the state estimation problem ---------------------
if strcmp('complex', sesettings.domain)
    results = run_cgn_sse(powsys, meas, sesettings);
elseif strcmp('real', sesettings.domain)
    if strcmp('sgn_sse', sesettings.method)
        results = run_gn_sse(powsys, meas, sesettings);
    elseif strcmp('wls_tse', sesettings.method)
        results = run_wls_tse_pmu(powsys, meas, sesettings);
    else
        return
    end
else
    return
end
% -------------------------------------------------------------------------

% -------------------------- Show results ---------------------------------
results_sse(sesettings, measurements, results)
% -------------------------------------------------------------------------
end

