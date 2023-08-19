function walk = randomwalk(samples)
walk = zeros(samples, 1);
for i = 2:samples
    walk(i) = walk(i - 1) + sign(randn);
end
end

