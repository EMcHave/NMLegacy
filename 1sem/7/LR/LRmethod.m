function [e1, V] = LRmethod(A)
Ao = A;
EPS = [];
K = [];
n = length(A);
V= [];
for t = 3:12
    eps = 10^-t;
    k = 0;
    A = Ao;
    while max(tril(A, -1), [], 'all')>eps
        [L,R] = LUfun(A);
        A = R*L;
        e1 = diag(A);
        k = k+1;
    end
    EPS = [EPS, eps];
    K = [K,k];
end

e1 = diag(e1);

for t = 1:n
    y = RevIt(Ao, e1(t,t)+0.1);
    v=y;
    V = [V, v];
end

figure
semilogx(EPS, K, '-r')
hold on
grid on
title('Зависимость кол-ва итераций от точности')
xlabel('Eps')
ylabel('Кол-во итераций')
end