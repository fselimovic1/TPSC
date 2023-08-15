function agcd = arraygcd(array)
agcd = array(1);
for i = 2:numel(array)
    agcd = gcd(agcd, array(i));
end
end

