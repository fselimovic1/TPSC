function [ col, elem ] = cJ_current_flow_phasor(nLine, data)
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
col =  [fromBus, toBus];
elem = [iiBranch, ijBranch];
end

