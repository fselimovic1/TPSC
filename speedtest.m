tic 
n = 10000;
a1 = randi(n, n);
a2 = randi(n, n);
a = [ a1 a2];
toc
tic
b = randi(n, n, 2 * n);
toc