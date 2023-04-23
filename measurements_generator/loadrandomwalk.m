function load = loadrandomwalk(samples, walk)
load = ones(1, samples);
for i = 2:samples
    load(i) = load(i - 1) + sign(randn) * walk;
end
end

