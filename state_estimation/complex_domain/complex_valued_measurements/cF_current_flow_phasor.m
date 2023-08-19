function h = cF_current_flow_phasor(nLine, branchi, branchj, ybus, state)
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
h = iiBranch * state(fromBus) + ijBranch * state(toBus);
end

