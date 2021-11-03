function Jd = JacobiDet(A, X0, eps)
n = length(A);
B = A*X0.';
D = zeros(n);
for i=1:n
    D(i,i) = A(i,i);
end

k = 0;
K = [];
Det = [];
E = [];

for i = 1:8
C = (eye(n)-D^-1*A);
g = D^-1*B;
q= norm(C, inf);
Xp = g;
Xn = C*Xp+g;
    while q/(1-q)*norm(Xn-Xp, inf)>eps
        t = Xn;
        Xp = Xn;
        Xn = C*t+g;
        k = k+1;
    end
K = [K, k];
k = 0;
Det = [Det, det(A)];
A(1,:)=2*A(1,:);
B(1) = 2*B(1);
D(1,1) = A(1,1);
end
figure 
semilogy(Det, K)
grid on
xlabel('Определитель')
ylabel('Количество итераций')
title('Зависимость кол-ва итераций от определителя матрицы')
Jd = K;
end