function agcd = arraygcd(array)
% adjuste to values in ms
array = array * 1000;
agcd = array(1);
for i = 2:numel(array)
    agcd = gcd(agcd, array(i));
end
agcd = agcd / 1000;
end

