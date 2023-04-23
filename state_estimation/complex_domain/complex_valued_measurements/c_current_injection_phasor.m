function [h, Hc] = c_current_injection_phasor(fromBus, data, state)
h = sparse([]);
if nargin == 3
    % measurements function
    h = sum(transpose(data.powerSystemAC.nodalMatrix(fromBus, :)) .* state(1:data.nBuses));
end

% Measurements jacobian matrix
Hc = transpose(data.powerSystemAC.nodalMatrix(fromBus, :));
end

