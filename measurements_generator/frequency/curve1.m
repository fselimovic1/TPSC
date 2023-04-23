function [fQ] = curve1(fM, fn, fMin, tS, maxT, dynT)

% quadraric function
t = [0, 1.58*maxT/2, maxT]; 
fQ = [fn, fMin, fMin + 0.15 * (fn - fMin)];
n = 2;
dT = 1/fM;
c = polyfit(t, fQ, n);
[t, fQ] = poly_fun(c, n, dT, max(t));

% superimpose dynamics
T = dynT/5;
kS = 0.3 * (fn - fMin);
kE = -(4 - dynT);
fQ = fQ - kS .* sin(2 * pi ./T .* t) .* exp(kE .* t);
fQ = [repmat(fn, 1, tS/dT), fQ];

% Smoothing the function
fc = 100; % cutoff frequency
fs = 1000; % sample rate
% Filter to achieve the smootheness
[b, a] = butter(7, fc/(fs/2), 'low'); % 6th order low-pass Butterworth filter
% Apply the filter
fQ = filter(b, a, fQ);
fQ(1:0.8 * fM) = fn;
end
