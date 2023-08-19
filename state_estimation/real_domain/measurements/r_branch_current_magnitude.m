function [ h, col, val ] = r_branch_current_magnitude(lineNo, data, state)
if lineNo > 0
    fromBus = data.branch(lineNo, 1);
    toBus = data.branch(lineNo, 2);
    iiBranch = data.powerSystemAC.nodalFromFrom(lineNo);
    ijBranch = data.powerSystemAC.nodalFromTo(lineNo);
else
    lineNo = -lineNo;
    fromBus = data.branch(lineNo, 2);
    toBus = data.branch(lineNo, 1);
    iiBranch = data.powerSystemAC.nodalToTo(lineNo);
    ijBranch = data.powerSystemAC.nodalToFrom(lineNo);
end
Vi = states(data.nBuses + fromBus);
Vj = states(data.nBuses + toBus);
thetaij = states(fromBus) - states(toBus);
hij = sqrt(data.powerSystemAC.A(lineNo) * states(data.nBuses + fromBus)^2 ...
     + data.powerSystemAC.B(lineNo) * states(data.nBuses + toBus)^2 - 2 * ...
     states(data.nBuses + toBus));
end

