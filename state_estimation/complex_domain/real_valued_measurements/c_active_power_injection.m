function [h, col, val] = c_active_power_injection(bus, ybus, state)
[~, col, nonZeros] = find(ybus.nodalMatrix(bus, :));
nBuses = size(ybus.nodalMatrix, 1);
% measurement function
h = state(bus) * (sum(conj(nonZeros) .* state(nBuses + col).'));
h = (h + conj(h)) / 2;

% measurements jacobian matrix
val = 1 / 2 .*  nonZeros .* state(nBuses + bus); 
busNo = find(col == bus);
val(busNo) = val(busNo) + 1 / 2 * (sum(conj(nonZeros) .* state(nBuses + col).'));
end



