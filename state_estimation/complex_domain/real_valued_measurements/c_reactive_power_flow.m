function [h, col, val] = c_reactive_power_flow(nLine, data, state)
nBuses = data.nBuses;
if nLine > 0
    fromBus = data.branch(nLine, 1);
    toBus = data.branch(nLine, 2);
    iiBranch = data.powerSystemAC.nodalFromFrom(nLine);
    ijBranch = data.powerSystemAC.nodalFromTo(nLine);
else
    fromBus = data.branch(-nLine, 2);
    toBus = data.branch(-nLine, 1);
    iiBranch = data.powerSystemAC.nodalToTo(-nLine);
    ijBranch = data.powerSystemAC.nodalToFrom(-nLine);
end


% measurement function
h = state(fromBus) * (state(fromBus + nBuses) * conj(iiBranch) + state(toBus + nBuses) * conj(ijBranch)); 
h = 1i * (conj(h) - h) / 2;

% measurements jacobian matrix
col = [ fromBus, toBus ];
val = 1i / 2 .* [ state(nBuses + fromBus) * iiBranch - state(nBuses + fromBus) * ...
        conj(iiBranch)  - state(nBuses + toBus) * conj(ijBranch), ...
        state(nBuses + fromBus) * ijBranch ];
end
