function [h, col, val] = c_active_power_injection(bus, data, state)
[~, col, nonZeros] = find(data.powerSystemAC.nodalMatrix(bus, :));
nBuses = data.nBuses;
% measurement function
h = state(bus) * (sum(conj(nonZeros) .* state(nBuses + col).'));
h = (h + conj(h)) / 2;

% measurements jacobian matrix
val = 1 / 2 .*  nonZeros .* state(nBuses + bus); 
busNo = find(col == bus);
val(busNo) = val(busNo) + 1 / 2 * (sum(conj(nonZeros) .* state(nBuses + col).'));
end



