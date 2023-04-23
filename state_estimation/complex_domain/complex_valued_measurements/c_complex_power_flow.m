function [h, Hc, Hconj] = c_complex_power_flow(nLine, data, state)
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

% measurements jacobian matrix
Hc = sparse(1, fromBus, state(fromBus + nBuses) * conj(iiBranch) + state(toBus + nBuses) * conj(ijBranch), 1, nBuses);
Hconj = sparse([1, 1], [fromBus, toBus], [state(fromBus) * conj(iiBranch), state(fromBus) * conj(ijBranch)], 1, nBuses);
end

