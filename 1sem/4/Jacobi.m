function Jacobi = Jacobi(A, X0, eps)
n = length(A);
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
Err = [];


while q/(1-q)*norm(Xn-Xp, inf)>eps
    t = Xn;
    Xp = Xn;
    Xn = C*t+g;
    k = k+1;
    K = [K, k];
    Err = [Err, norm(X0.'-Xn)];
end
figure 
semilogy(K, Err)
grid on
xlabel('K')
ylabel('Ошибка')
title('Зависимость ошибки от номера итерации')
Jacobi = X;
end

