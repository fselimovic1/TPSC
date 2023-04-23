function load = sudden_change(dS, start, range, load)
step = dS/range;
acc = 0;
for i = start:start + range
    acc = acc + step;
    load(i) = load(i) + acc;
end
load(i + 1:end) = load(i + 1:end) + acc;
end

