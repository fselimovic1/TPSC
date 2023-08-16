function [ col, elem ] = cJ_current_flow_phasor(nLine, branchi, branchj, ybus)
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
col =  [ fromBus, toBus ];
elem = [ iiBranch, ijBranch ];
end

