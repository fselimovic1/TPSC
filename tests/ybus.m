clc
clear 
 

ntrails = 50;

casename = 'case13659pegase';

data = loadcase(casename);

% number of:
num.bus = size(data.bus, 1);
num.branch = size(data.branch, 1);
num.gen = size(data.generator, 1);

% read the power system data
bus = 1:num.bus;
busi = data.bus(:, 1);
branchi = data.branch(:, 1);
branchj = data.branch(:, 2);
resistance = data.branch(:, 3);
reactance = data.branch(:, 4);
charging = data.branch(:, 5);
transturns = data.branch(:, 9);
transphase = data.branch(:, 10);
branchstatus = data.branch(:, 11) == 1;
bustype = data.bus(:, 2);
Pload = data.bus(:, 3);
Qload = data.bus(:, 4);
Gshunt = data.bus(:, 5);
Bshunt = data.bus(:, 6);
Vmi = data.bus(:, 8);
Vai = data.bus(:, 9);
genbuses = data.generator(:, 1);
Pgeni = data.generator(:, 2);
Qgeni = data.generator(:, 3);
% Qgmax, Qgmin
Vgen = data.generator(:, 6);

num.islack= find(bustype == 3);

times = zeros(ntrails, 1);
% ps1 = admittance_matrix(num, bus, branchi, branchj, resistance, reactance, charging, transturns, transphase, branchstatus, Gshunt, Bshunt);
% ps2 = admittance_matrix2(num, bus, branchi, branchj, resistance, reactance, charging, transturns, transphase, branchstatus, Gshunt, Bshunt);
% sum(sum(abs(ps1.nodalMatrix - ps2.nodalMatrix) > 1e-6))
tic
for i = 1:ntrails
    powerSystem = admittance_matrix2(num, bus, branchi, branchj, resistance, reactance, charging, transturns, transphase, branchstatus, Gshunt, Bshunt);
end
times(i) = toc;
fprintf("Time required to compute ybus matrix for the %s is %.2f [ms].\n", casename, mean(times) * 1000);
