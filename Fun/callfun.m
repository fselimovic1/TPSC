clc
clear 

fSubType = 2; % / 2
if fSubType == 1
    tS = 1;
    tSE = 5;
    fn = 50;
    fMin = 49;
    fM = 100;
    tDyn = 0.2 * tSE;
freq_profile = curve1(fM, fn, fMin, tS, tSE - tS, tDyn);
elseif fSubType == 2
    tS = 1;
    tE = 0;
    fn = 50;
    tSE = 5;
    fMin = 49;
    fM = 100;
    freq_profile = curve2(fM, tS, tE, tSE - tS - tE, fn, fMin);
end

% PLOT results
t = 0:1/fM:tSE - 1/fM;
plot(t, freq_profile)
xlabel('t [s]')
ylabel('f [Hz]')