function [fC] = curve2(fM, tS, tE, tCon, fn, fMin)
% INPUTS:
% fM - maximal PMU reporting frequency;
% fn - nominal frequency
% fMin - frequency minimum on curve
% tS - steady-state time (start)
% tE - steady-state time (end)
% tCon - time of dynamics


dT = 1/fM;
% data from graph
tD = [ 0, 5, 15, 20, 30,  40] .* (tCon/40); 
fD = [ 0, -1, -0.73, -0.75, -0.4, -0.1] .* (fn - fMin);
t = 0:dT:tCon;
f = interp1(tD, fD, t, 'spline');
n = 20;
c = polyfit(t, f, n);
[~, f] = poly_fun(c, n, dT, tCon);
fC = [repmat(fn, 1, (tS + tCon)/dT + 1), repmat(fn + f(end), 1, tE/dT)];
fC(tS/dT + 1:(tS + tCon)/dT + 1 ) = fC(tS/dT + 1:(tS + tCon)/dT + 1) + f;
end

