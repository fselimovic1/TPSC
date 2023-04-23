function [t, fun] = poly_fun(c, n, dT, T)
fun = 0;
t = 0:dT:T - dT; 
for i = 1 : n + 1
    fun = fun + c(i) .* t .^(n - i + 1);
end
end

