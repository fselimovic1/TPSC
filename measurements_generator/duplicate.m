function boolean = duplicate(k, l, set)
boolean = 0;
for i = 1:size(set, 1)
    if set(i, 2) == k && set(i, 3) == l
        boolean = 1;
    end
end
end

