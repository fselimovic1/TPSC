function [ results, data ] = runpf(casename, pfsettings, varargin)
% load power system case
data = loadcase(casename);

% number of elements in power system
num.bus = size(data.bus, 1);
num.branch = size(data.branch, 1);
num.gen = size(data.generator, 1);
num.islack = find(data.bus(:, 2) == 3);

% load power system data into arrays 
powersystem.bus = 1:num.bus;
powersystem.busi = data.bus(:, 1);
powersystem.branchi = data.branch(:, 1);
powersystem.branchj = data.branch(:, 2);
powersystem.resistance = data.branch(:, 3);
powersystem.reactance = data.branch(:, 4);
powersystem.charging = data.branch(:, 5);
powersystem.transturns = data.branch(:, 9);
powersystem.transphase = data.branch(:, 10);
powersystem.branchstatus = data.branch(:, 11) == 1;
powersystem.bustype = data.bus(:, 2);
powersystem.Pload = data.bus(:, 3);
powersystem.Qload = data.bus(:, 4);
powersystem.Gshunt = data.bus(:, 5);
powersystem.Bshunt = data.bus(:, 6);
powersystem.Vmi = data.bus(:, 8);
powersystem.Vai = data.bus(:, 9);
powersystem.genbuses = data.generator(:, 1);
powersystem.Pgeni = data.generator(:, 2);
powersystem.Qgeni = data.generator(:, 3);
% Qgmax, Qgmin
powersystem.Vgen = data.generator(:, 6);
% bus types
powersystem.isPQ = powersystem.bustype == 1;
powersystem.isPV = powersystem.bustype == 2;

% renumber data if needed
toRenumber = any(data.bus(:, 1) ~= powersystem.bus);
if toRenumber
    powersystem.genbuses = renumbering(powersystem.genbuses, powersystem.busi, powersystem.bus);
    powersystem.branchi = renumbering(powersystem.branchi, powersystem.busi, powersystem.bus);
    powersystem.branchj = renumbering(powersystem.branchj, powersystem.busi, powersystem.bus);
end

if nargin == 3
    dynsettings = varargin{1};
    % frequency - realation with model parameters
    powersystem.reactance = powersystem.reactance * dynsettings.f/data.fn;
    powersystem.charging = powersystem.charging * data.fn/dynsettings.f;
    powersystem.Bshunt = powersystem.Bshunt * data.fn/dynsettings.f;
    
    % load dynamics
    if dynsettings.loadNo
        if dynsettings.loadNo == -1
            powersystem.Pload = powersystem.Pload .* dynsettings.load;
            powersystem.Qload = powersystem.Qload .* dynsettings.load;
        else
            powersystem.Pload(dynsettings.loadNo) = powersystem.Pload(dynsettings.loadNo) + ...
                                   dynsettings.load(1); 
            powersystem.Qload(dynsettings.loadNo) = powersystem.Qload(dynsettings.loadNo) + ...
                                   dynsettings.load(2);                   
        end
    end
end
tic
% power system admittance matrix
ybus = admittance_matrix(num, powersystem.bus, powersystem.branchi, powersystem.branchj, ...
    powersystem.resistance, powersystem.reactance, powersystem.charging,...
    powersystem.transturns, powersystem.transphase, powersystem.branchstatus, ...
    powersystem.Gshunt, powersystem.Bshunt);

% Solve the power flows analysis problem
if strcmp('complex', pfsettings.domain)
	[ Vc, iter, converged, method ] = run_cnr_pf(pfsettings, powersystem, ybus, num);
end
algtime = toc;
% Post processing of results
if converged
    if nargin == 3
        Vc = abs(Vc(powersystem.bus)) .* exp(1i .* (angle(Vc(powersystem.bus)) + dynsettings.thetaslack));
    end
    results = postprocess_acpf(ybus, powersystem.branchi, powersystem.branchj, Vc);
    results.Pload = powersystem.Pload;
    results.Qload = powersystem.Qload;
    results.Pgen(powersystem.genbuses) = powersystem.Pgeni;
    results.Qgen(powersystem.genbuses) = results.Qi(powersystem.genbuses)...
                                         + results.Qload(powersystem.genbuses);
end
results.converged = converged;
results.iter = iter;
results.algtime = algtime;
results.method = method;
results.sys = casename;
if pfsettings.info
    results_pf(pfsettings, results, num, powersystem.bus, powersystem.branchi, powersystem.branchj, data.baseMVA);
end
end
