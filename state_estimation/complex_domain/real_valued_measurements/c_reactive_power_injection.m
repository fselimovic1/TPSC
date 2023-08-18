function [h, col, val] = c_reactive_power_injection(bus, ybus, state)
[~, col, nonZeros] = find(ybus.nodalMatrix(bus, :));
nBuses = size(ybus.nodalMatrix, 1);
% measurement function
h = state(bus) * (sum(conj(nonZeros) .* state(nBuses + col).'));
h = 1i * (conj(h) - h) / 2;

% measurements jacobian matrix
val =  nonZeros .* state(nBuses + bus); 
busNo = find(col == bus);
val(busNo) = val(busNo) - sum(conj(nonZeros) .* state(nBuses + col).');
val = (1i/2) .* val; 
end
