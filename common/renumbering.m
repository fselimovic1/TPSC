function array = renumbering(array, original, new)
pairs = containers.Map(original, new);
for i = 1:numel(array)
    array(i) = pairs(array(i));
end
end

