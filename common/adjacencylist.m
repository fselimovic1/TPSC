function adj = adjacencylist(nBuses, branchi, branchj)
    adj = cell(nBuses, 1);
    for i = 1:nBuses
        adj{i} = unique([branchj(branchi == i); branchi(branchj == i)]);
    end
end


