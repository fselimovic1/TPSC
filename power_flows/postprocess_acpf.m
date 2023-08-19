function [ results ] = postprocess_acpf(ybus, branchi, branchj, v)
nBranches = numel(ybus.admittance);
% bus voltage
results.Vm = abs(v);
results.Va = angle(v);

% injected current
Ii = ybus.nodalMatrix * v;
results.Iim = abs(Ii);
results.Iia = angle(Ii);

% injected powers
Si = v .* conj(Ii);
results.Pi = real(Si);
results.Qi = imag(Si);
results.Pgen = zeros(numel(Si), 1);
results.Qgen = results.Pgen;
results.Pload = results.Pgen;
results.Qload = results.Pgen;

% current flows
results.Iijm = zeros(nBranches, 1);
results.Ijim = zeros(nBranches, 1);
results.Iija = zeros(nBranches, 1);
results.Ijia = zeros(nBranches, 1);
results.Pij  = zeros(nBranches, 1);
results.Pji  = zeros(nBranches, 1);
results.Qij  = zeros(nBranches, 1);
results.Qji  = zeros(nBranches, 1);
results.Ploss = zeros(nBranches, 1);
results.Qloss = zeros(nBranches, 1);

for i = 1:nBranches
    Iij = [ ybus.nodalFromFrom(i), ybus.nodalFromTo(i)
           ybus.nodalToFrom(i), ybus.nodalToTo(i)] *...
           [ v(branchi(i)); v(branchj(i))];
    Sij = v(branchi(i)) * conj(Iij(1));
    Sji = v(branchj(i)) * conj(Iij(2));
    results.Iijm(i) = abs(Iij(1));
    results.Ijim(i) = abs(Iij(2));
    results.Iija(i) = angle(Iij(1));
    results.Ijia(i) = angle(Iij(2));
    results.Pij(i)  = real(Sij);
    results.Pji(i)  = real(Sji);
    results.Qij(i)  = imag(Sij);
    results.Qji(i)  = imag(Sji);
    
    % losses on the branch
    Ib = ybus.admittance(i) * (v(branchi(i)) / ...
        ybus.transformerRatio(i) - v(branchj(i)));
    results.Ploss(i) = abs(Ib)^2 * real(1 / ybus.admittance(i));
    results.Qloss(i) = abs(Ib)^2 * imag(1 / ybus.admittance(i));
end
end