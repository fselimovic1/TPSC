function [h, col, val] = c_branch_current_magnitude(lineNo, data, state)
if lineNo > 0
    fromBus = data.branch(lineNo, 1);
    toBus = data.branch(lineNo, 2);
    iiBranch = data.powerSystemAC.nodalFromFrom(lineNo);
    ijBranch = data.powerSystemAC.nodalFromTo(lineNo);
else
    fromBus = data.branch(-lineNo, 2);
    toBus = data.branch(-lineNo, 1);
    iiBranch = data.powerSystemAC.nodalToTo(-lineNo);
    ijBranch = data.powerSystemAC.nodalToFrom(-lineNo);
end
ijPh = iiBranch * state(fromBus) + ijBranch * state(toBus);
h = sqrt(ijPh * conj(ijPh));
col = [ fromBus, toBus ];
ijAngle =  angle(ijPh);
val = [iiBranch / 2 * exp(-1i * ijAngle), ijBranch / 2 * exp(-1i * ijAngle)];
end

