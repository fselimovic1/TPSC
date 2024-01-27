function walk = randomwalk(samples)
choices = [ -1, 0, 1 ];
walk = zeros(samples, 1);
for i = 2:samples
    walk(i) = walk(i - 1) + choices(randi(3));
end
end

