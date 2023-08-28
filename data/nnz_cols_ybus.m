function cols = nnz_cols_ybus(nBuses, branchi, branchj)
    cols = cell(nBuses, 1);
    for i = 1:nBuses
        cols{i} = sort([i; unique([branchj(branchi == i); branchi(branchj == i)])]);
    end
end
