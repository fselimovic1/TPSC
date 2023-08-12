function data = renumbering(data)
for i = 1:size(data.bus, 1)
    if i ~= data.bus(i, 1)
        data.branch(data.branch(:, 1) == data.bus(i, 1), 1) = i;
        data.branch(data.branch(:, 2) == data.bus(i, 1), 2) = i;
        data.gen(data.gen(:, 1) == data.bus(i, 1), 1) = i;
        data.bus(i, 1) = i;
    end
end
end

