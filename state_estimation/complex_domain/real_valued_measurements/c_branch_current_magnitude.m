function [h, col, val] = c_branch_current_magnitude(nLine, branchi, branchj, ybus, state)
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
ijPh = iiBranch * state(fromBus) + ijBranch * state(toBus);
h = sqrt(ijPh * conj(ijPh));
col = [ fromBus, toBus ];
ijAngle =  angle(ijPh);
val = [iiBranch / 2 * exp(-1i * ijAngle), ijBranch / 2 * exp(-1i * ijAngle)];
end

