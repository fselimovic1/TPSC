function [ results, data ] = runpf(casename, pfsettings, varargin)
% ----------------------- Load Power System -------------------------------
data = loadcase(casename);
%--------------------------------------------------------------------------

%--------------------- Extract Useful Informations ------------------------
powsys = preprocess_ps(data, 'pf');
%--------------------------------------------------------------------------

% ---------------------- Dynamics conditions (runmg) ----------------------
if nargin == 3
    dynsettings = varargin{1};
    %------------- Frequency - realation with model parameters ------------
    powsys.branch.reactance = powsys.branch.reactance * dynsettings.f/data.fn;
    powsys.branch.charging = powsys.branch.charging * data.fn/dynsettings.f;
    powsys.bus.Bshunt = powsys.bus.Bshunt * data.fn/dynsettings.f;
    % ---------------------------------------------------------------------
    
    % ------------------------- Load dynamics -----------------------------
    if dynsettings.loadNo
        if dynsettings.loadNo == -1
            powsys.bus.Pload = powsys.bus.Pload .* dynsettings.load;
            powsys.bus.Qload = powsys.bus.Qload .* dynsettings.load;
        else
            powsys.bus.Pload(dynsettings.loadNo) = powsys.bus.Pload(dynsettings.loadNo) + ...
                                   dynsettings.load(1); 
            powsys.bus.Qload(dynsettings.loadNo) = powsys.bus.Qload(dynsettings.loadNo) + ...
                                   dynsettings.load(2);                   
        end
    end
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

% --------------------- Calculate Y matrix --------------------------------
powsys = admittance_matrix(powsys);
% -------------------------------------------------------------------------

% ----------------- Solve the power flows analysis problem ----------------
tic
if strcmp('complex', pfsettings.domain)
	[ Vc, iter, converged, method ] = run_cnr_pf(powsys, pfsettings);
end
algtime = toc;
% -------------------------------------------------------------------------

% --------------------- Post processing of the results --------------------
if converged
    if nargin == 3
        Vc = abs(Vc) .* exp(1i .* (angle(Vc) + dynsettings.thetaslack));
    end
    results = postprocess_acpf(powsys, Vc);
end
results.converged = converged;
results.iter = iter;
results.algtime = algtime;
results.method = method;
results.sys = casename;
% -------------------------------------------------------------------------

if pfsettings.info
    results_pf(results, powsys, pfsettings, data.baseMVA);
end
end
