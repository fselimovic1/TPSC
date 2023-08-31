function [ results, measurements, data ] = runsse(casename, vrs, sesettings)
% ----------------- Load Power System Data and Measurements ---------------
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');
% -------------------------------------------------------------------------

%----------------- Extract Useful Informations (Power System) -------------
powsys = preprocess_ps(data, 'se');
%--------------------------------------------------------------------------

% --------------------- Calculate Y matrix --------------------------------
powsys = admittance_matrix(powsys);
% -------------------------------------------------------------------------

%----------------- Extract Useful Informations (Measurements) -------------
meas = preprocess_meas(data, measurements);
%--------------------------------------------------------------------------

% ---------------- Solve the state estimation problem ---------------------
tic
if strcmp('complex', sesettings.domain)
    if strcmp('cgn_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, sesettings);
    elseif strcmp('cls_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_cls_sse(powsys, meas, sesettings);
    end
else
    return
end
algtime = toc;
% -------------------------------------------------------------------------

% -------------------------- Show results ---------------------------------
results.voltage = Vc;
results.converged = converged;
results.iter = iter;
results.algtime = algtime;
results.info = info;
results.sys = casename;

results_sse(results, powsys, meas, sesettings)
% -------------------------------------------------------------------------
end

