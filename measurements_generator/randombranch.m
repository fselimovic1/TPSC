function [busidx, idx] = randombranch(branch, nBranches, nM)
randidx = randperm(2 * nBranches);
randidx(randidx > nBranches) = randidx(randidx > nBranches) - 2 * nBranches - 1;
idx = randidx(1:nM)';
busidx = zeros(nM, 1);
busidx(idx < 0) = branch(abs(idx(idx < 0)), 2);
busidx(idx > 0) = branch(idx(idx > 0), 1);

end

