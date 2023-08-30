function [ results ] = postprocess_acpf(powsys, Vc)
% --------------------------- Bus voltage ---------------------------------
results.Vm = abs(Vc);
results.Va = angle(Vc);
%--------------------------------------------------------------------------

% ------------------------- Injected current ------------------------------
Ii = powsys.ybus.y * Vc;
results.Iim = abs(Ii);
results.Iia = angle(Ii);
% -------------------------------------------------------------------------

% ------------------------- Injected powers -------------------------------
Si = Vc .* conj(Ii);
results.Pi = real(Si);
results.Qi = imag(Si);
results.Pload = powsys.bus.Pload;
results.Qload = powsys.bus.Qload;
results.Pgen = powsys.bus.Pgen;
results.Qgen = zeros(powsys.num.bus, 1);
results.Qgen(powsys.gen.bus) = results.Qi(powsys.gen.bus)...
                                         + results.Qload(powsys.gen.bus);
% -------------------------------------------------------------------------

% ------------------- Branch Current and Power Flows ----------------------
Iij = [ powsys.ybus.fromfrom, powsys.ybus.fromto ] * [ Vc(powsys.branch.i), Vc(powsys.branch.j)].';
Iij = diag(Iij);
Iji = [ powsys.ybus.tofrom, powsys.ybus.toto ] * [ Vc(powsys.branch.i), Vc(powsys.branch.j)].';
Iji = diag(Iji);
Sij = Vc(powsys.branch.i) .* conj(Iij);
Sji = Vc(powsys.branch.j) .* conj(Iji);
results.Iijm = abs(Iij);
results.Ijim = abs(Iji);
results.Iija = angle(Iij);
results.Ijia = angle(Iji);
results.Pij  = real(Sij);
results.Pji  = real(Sji);
results.Qij  = imag(Sij);
results.Qji  = imag(Sji);
% -------------------------------------------------------------------------

% ---------------------------- Losses on branches -------------------------
Ib = powsys.ybus.admittance .* (Vc(powsys.branch.i) ./ powsys.ybus.transratio - Vc(powsys.branch.j));
results.Ploss = abs(Ib).^2 .* real(1 ./ powsys.ybus.admittance);
results.Qloss = abs(Ib).^2 .* imag(1 ./ powsys.ybus.admittance);
% -------------------------------------------------------------------------
end