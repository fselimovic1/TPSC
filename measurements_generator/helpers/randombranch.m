function [busidx, idx] = randombranch(branchi, branchj, nBranches, nM)
randidx = randperm(2 * nBranches);
randidx(randidx > nBranches) = randidx(randidx > nBranches) - 2 * nBranches - 1;
idx = randidx(1:nM)';
busidx = zeros(nM, 1);
busidx(idx < 0) = branchi(abs(idx(idx < 0)));
busidx(idx > 0) = branchj(idx(idx > 0));
end

