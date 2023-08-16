function adj = adjacencylist(nBuses, branchi, branchj)
adj = cell(nBuses, 1);
for i = 1:nBuses
    adj{i} = branchj(branchi == i);
    adj{i} = [adj{i}; branchi(branchj == i)];
    adj{i} = unique(adj{i});  
end
end

