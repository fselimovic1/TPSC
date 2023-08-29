function [ results, data ] = runpf(casename, pfsettings, varargin)
% ----------------------- Load Power System -------------------------------
data = loadcase(casename);
%--------------------------------------------------------------------------
tic
%--------------------- Extract Useful Informations ------------------------
powsys = preprocess(data, 'pf');
%--------------------------------------------------------------------------

% ---------------------- Dynamics conditions (runmg) ----------------------
if nargin == 3
    dynsettings = varargin{1};
    % frequency - realation with model parameters
    powsys.branch.reactance = powsys.branch.reactance * dynsettings.f/data.fn;
    powsys.branch.charging = powsys.charging * data.fn/dynsettings.f;
    powsys.bus.Bshunt = powsys.bus.Bshunt * data.fn/dynsettings.f;
    
    % load dynamics
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
end
% -------------------------------------------------------------------------
% --------------------- Calculate Y matrix --------------------------------
powsys = admittance_matrix(powsys);
% -------------------------------------------------------------------------
toc,tic
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
        Vc = abs(Vc(powsys.bus)) .* exp(1i .* (angle(Vc(powsys.bus)) + dynsettings.thetaslack));
    end
    results = postprocess_acpf(ybus, powsys.branchi, powsys.branchj, Vc);
    results.Pload = powsys.Pload;
    results.Qload = powsys.Qload;
    results.Pgen(powsys.genbuses) = powsys.Pgeni;
    results.Qgen(powsys.genbuses) = results.Qi(powsys.genbuses)...
                                         + results.Qload(powsys.genbuses);
end
results.converged = converged;
results.iter = iter;
results.algtime = algtime;
results.method = method;
results.sys = casename;
% -------------------------------------------------------------------------

if pfsettings.info
    results_pf(pfsettings, results, num, powsys.bus, powsys.branchi, powsys.branchj, data.baseMVA);
end
end
