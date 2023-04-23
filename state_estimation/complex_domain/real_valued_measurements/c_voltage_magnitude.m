function [h, H] = c_voltage_magnitude(bus, nBuses, state)
% measurement function
h = sqrt(state(bus) * state(bus + nBuses));

% measurements jacobian matrix
H = sparse(1, bus, 0.5 * exp(-1i * angle(state(bus))), 1, nBuses);
end

