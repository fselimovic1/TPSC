function [h, Hc] = c_voltage_phasor(fromBus, nBuses, x)
if nargin == 3
    h = x(fromBus);
end
Hc = sparse(1, nBuses);
Hc(fromBus) = 1;
end

