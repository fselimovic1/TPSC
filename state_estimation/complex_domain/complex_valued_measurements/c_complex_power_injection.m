function [h, Hc, Hconj] = c_complex_power_injection(bus, data, state)
nBuses = data.nBuses;
% measurement function
h = state(bus) * sum(conj(data.powerSystemAC.nodalMatrix(bus, :)) .* state(nBuses + 1:2 * data.nBuses).');

% measurements jacobian matrix
Hc = sparse(1, bus, sum(conj(data.powerSystemAC.nodalMatrix(bus, :)) .* state(nBuses + 1:2 * nBuses).'), 1, nBuses);
[~, cols, vals] = find(data.powerSystemAC.nodalMatrix(bus, :));
nNonZero = numel(vals);
Hconj = sparse(ones(1, nNonZero), cols, conj(vals) .* state(bus), 1, nBuses);
end

