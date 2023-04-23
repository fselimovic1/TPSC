function [ load ] = load_profile(fM, tSE)
x = 0:tSE/96:tSE;
x = x(1:end - 1);
t = 0:1/fM:tSE;
t = t(1:end - 1);
load = interp1(x, xlsread('load_data.xlsx', 'Sheet1', 'B2:B97'), t, 'pchip');
load = load./max(load);
end

