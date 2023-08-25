function [ results, measurements, data ] = runrtse(casename, vrs, rtsesettings)
% --------------------- Load PS Data and Measurements ---------------------
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');



end

