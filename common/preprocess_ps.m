function powsys = preprocess_ps(data, type)
% -------------------- Number of elements in power sytem ------------------
branchstatus = logical(data.branch(:, 11));
genstatus = logical(data.generator(:, 8));
powsys.num.bus = size(data.bus, 1);
powsys.num.branch = sum(branchstatus);
powsys.num.gen = sum(genstatus);
powsys.num.islack = find(data.bus(:, 2) == 3);
% -------------------------------------------------------------------------
powsys.fn = data.fn;
% -------------------------- Bus data -------------------------------------
powsys.bus.busnew = (1:powsys.num.bus)';
powsys.bus.busorg = data.bus(:, 1);
powsys.bus.Vmi = data.bus(:, 8);
powsys.bus.Vai = data.bus(:, 9);
powsys.bus.Gshunt = data.bus(:, 5);
powsys.bus.Bshunt = data.bus(:, 6);
if strcmp(type, 'pf')
    powsys.bus.type = data.bus(:, 2);
    powsys.bus.Pload = data.bus(:, 3);
    powsys.bus.Qload = data.bus(:, 4);
    powsys.bus.Vmin = data.bus(:, 13);
    powsys.bus.Vmax = data.bus(:, 12);
end
if strcmp(type, 'se')
    % ------------------------- ZI buses ----------------------------------
    powsys.bus.zi = find(data.bus(:, 2) == 1 & data.bus(:, 3) == 0 & data.bus(:, 4) == 0);
    powsys.num.zi = numel(powsys.bus.zi); 
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

% -------------------------- Branch data ----------------------------------
powsys.branch.i = data.branch(branchstatus, 1);
powsys.branch.j = data.branch(branchstatus, 2);
powsys.branch.resistance = data.branch(branchstatus, 3);
powsys.branch.reactance = data.branch(branchstatus, 4);
powsys.branch.charging = data.branch(branchstatus, 5);
powsys.branch.transturns = data.branch(branchstatus, 9);
powsys.branch.transphase = data.branch(branchstatus, 10);
% -------------------------------------------------------------------------

% ------------------------- Generator -------------------------------------
powsys.gen.bus = data.generator(genstatus, 1);
% -------------------------------------------------------------------------

% ---------------------- Renumber data if needed --------------------------
toRenumber = any(powsys.bus.busnew ~= powsys.bus.busorg);
if toRenumber
    powsys.gen.bus = renumbering(powsys.gen.bus, powsys.bus.busorg,  powsys.bus.busnew);
    powsys.branch.i = renumbering(powsys.branch.i, powsys.bus.busorg,  powsys.bus.busnew);
    powsys.branch.j = renumbering(powsys.branch.j, powsys.bus.busorg,  powsys.bus.busnew);
end
% -------------------------------------------------------------------------

% -------------------------- Generator -> Bus -----------------------------
A = sparse(powsys.gen.bus, 1:powsys.num.gen, 1, powsys.num.bus, powsys.num.gen);
if strcmp(type, 'pf')
    powsys.bus.Pgen = A * data.generator(genstatus, 2);
    powsys.bus.Qgen = A * data.generator(genstatus, 3);
    powsys.bus.Qgmin = A * data.generator(genstatus, 10);
    powsys.gen.Qgmax = A * data.generator(genstatus, 9);
    powsys.bus.Vmi(powsys.gen.bus) = data.generator(genstatus, 6);
    if ~all(genstatus)
        powsys.bus.type = 1;
        powsys.bus.type(sparse(powsys.gen.bus, 1, 1, ...
        powsys.num.bus, 1) .* data.bus(:, 2) == 2) = 2;
        powsys.bus.type(powsys.num.islack) = 3;
    end
    powsys.bus.pq = powsys.bus.busnew(powsys.bus.type == 1);
    powsys.bus.pv = powsys.bus.busnew(powsys.bus.type == 2);    
end
% -------------------------------------------------------------------------

if strcmp("se", type)
    % ----------------- PMU Measurements Devices --------------------------
    powsys.pmu.onbus = data.pmu(:, 1);
    powsys.pmu.rfreq = data.pmu(:, 7);
    powsys.pmu.msd = data.pmu(:, 3);
    powsys.pmu.asd = data.pmu(:, 4);
    powsys.pmu.fsd = data.pmu(:, 5);
    powsys.pmu.rfsd = data.pmu(:, 6);
    % ---------------------------------------------------------------------
    
    % -------------------- SCADA Measurement Devices ----------------------
    powsys.scada.type = data.pmu(:, 2);
    powsys.scada.onbus = data.pmu(:, 1);
    powsys.pmu.rfreq = data.pmu(:, 5);
    powsys.pmu.sd = data.pmu(:, 4);
    % ---------------------------------------------------------------------
end
end

