function J = JacobiEps(A, X0)
n = length(A);
eps0 = 10^-3;
B = A*X0.';
D = zeros(n);
for i=1:n
    D(i,i) = A(i,i);
end


C = (eye(n)-D^-1*A);
g = D^-1*B;
q = norm(C, inf);

Xp = g;
Xn = C*Xp+g;
k = 0;
K = [];
Eps = [];
delta = [];
eps = eps0;
for i = 3:13
    eps = 10^-i;
    while q/(1-q)*norm(Xn-Xp, inf)>eps
        t = Xn;
        Xp = Xn;
        Xn = C*t+g;
        k = k+1;
    end
    K = [K, k];
    Eps = [Eps, eps];
    delta = [delta, norm(X0.'-Xn, inf)];
end
figure 
semilogx(Eps, K)
grid on
xlabel('Eps')
ylabel('К')
title('Зависимость количества итераций от заданной точности')

figure 
loglog(Eps, delta)
grid on
xlabel('Eps')
ylabel('||X0-X||')
title('Зависимость погрешности от заданной точности')
J = delta;
end
