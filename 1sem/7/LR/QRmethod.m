function [e1, V] = QRmethod(A, H, q)
Ao = A;
Ho = H;
EPS = [];
K = [];
K1 = [];
n = length(A);
V= [];
for t = 3:10
    eps = 10^-t;
    k = 0;
    k1 = 0;
    A = Ao;
    H = Ho;
    while max(tril(A, -1), [], 'all')>eps
        [Q,R] = givens(A);
        A = Q.'*A*Q;
        e1 = diag(A);
        k = k+1;
    end
    while max(tril(H, -1), [], 'all')>eps
        [Q,R] = givensMod(H);
        H = Q.'*H*Q;
        eh1 = diag(H);
        k1 = k1+1;
    end
    EPS = [EPS, eps];
    K = [K,k];
    K1 = [K1, k1] ; 

end

e1 = diag(e1);

for t = 1:n
    y = RevIt(Ho, e1(t,t)+0.1);
    v=q*y;
    V = [V, v];
end

figure
semilogx(EPS, K, '-r')
hold on
grid on
semilogx(EPS, K1, '-b')
title('Зависимость кол-ва итераций от точности(хорошая отделимость)')
legend('Изначальная матрица', 'Форма Хессенберга')
xlabel('Eps')
ylabel('Кол-во итераций')
end