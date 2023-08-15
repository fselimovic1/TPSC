function [fQ] = curve1(fM, fn, fMin, tS, maxT, dynT)
% INPUTS:
% fM - maximal PMU reporting frequency;
% fn - nominal frequency
% fMin - frequency minimum on curve
% tS - steady-state time
% maxT - time of dynamics
% dynT - time of oscillations


% quadraric function
t = [0, 1.58*maxT/2, maxT]; 
fQ = [fn, fMin, fMin + 0.15 * (fn - fMin)];
n = 2;
dT = 1/fM;
c = polyfit(t, fQ, n);
[t, fQ] = poly_fun(c, n, dT, max(t));

% + dynamics
T = dynT/5;
kS = 0.3 * (fn - fMin);
kE = -(4 - dynT);
fQ = fQ - kS .* sin(2 * pi ./T .* t) .* exp(kE .* t);
fQ = [repmat(fn, 1, tS/dT), fQ];
end
