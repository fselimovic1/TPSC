function [ idxa, idxb ] = myintersect(array, aprop, bprop)
idxa = [];
idxb = [];
if ~isempty(aprop) && ~isempty(bprop)
    for i = 1:numel(array)
        [ ~, ib ] = ismember(array(i), array(i+1:end));
        if ismember(i, aprop) && ismember(i + ib, bprop)
            idxa = [ idxa; i ];
            idxb = [ idxb; i + ib ];
        end
    end
end
end

