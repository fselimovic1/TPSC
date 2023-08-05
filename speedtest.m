% n = 100000000;
% tic
% a1 = complex(zeros(n, 1));
% a2 = complex(zeros(n, 1));
% a3 = complex(zeros(n, 1));
% a4 = complex(zeros(n, 1));
% a5 = complex(zeros(n, 1));
% toc
% tic
% b1 = complex(zeros(n, 1));
% b2 = b1;
% b3 = b1;
% b4 = b1;
% b5 = b1;
% toc
% b1(20) = 20;
% fprintf("%d, %d", b1(20), b3(20))
n = 10000;
tic
f1(randi(50, n))
toc
tic
f2(randi(50, n))
toc
function f1(A)
[m, n] = size(A);
mi = randi(m);
ni = randi(n);
A(mi, ni) = randi(m);
fprintf("Num: %d\n", A(mi, ni))
end 
function f2(A)
[m, n] = size(A);
mi = randi(m);
ni = randi(n);
fprintf("Num: %d\n", A(mi, ni))
end 