function [ h, col, val ] = c_active_power_flow(nLine, branchi, branchj, ybus, state)
nBuses = size(ybus.nodalMatrix, 1);
if nLine > 0
    fromBus = branchi(nLine);
    toBus = branchj(nLine);
    iiBranch = ybus.nodalFromFrom(nLine);
    ijBranch = ybus.nodalFromTo(nLine);
else
    fromBus = branchj(-nLine);
    toBus = branchi(-nLine);
    iiBranch = ybus.nodalToTo(-nLine);
    ijBranch = ybus.nodalToFrom(-nLine);
end

% measurement function
h = state(fromBus) * (state(fromBus + nBuses) * conj(iiBranch) + state(toBus + nBuses) * conj(ijBranch)); 
h = (h + conj(h)) / 2;

% measurements jacobian matrix
col = [ fromBus, toBus ];
val = 1 / 2 .* [ state(nBuses + fromBus) * conj(iiBranch) + state(nBuses + toBus) * ...
        conj(ijBranch) + state(nBuses + fromBus) * iiBranch, ...
        state(nBuses + fromBus) * ijBranch ];
end

