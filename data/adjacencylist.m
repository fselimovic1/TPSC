function adj = adjacencylist(data)
adj = cell(data.nBuses, 1);
for i = 1:data.nBuses
    adj{i} = data.branch(data.branch(:, 1) == i, 2)';
    adj{i} = [adj{i},  data.branch(data.branch(:, 2) == i, 1)'];
    adj{i} = unique(adj{i});  
end
end

