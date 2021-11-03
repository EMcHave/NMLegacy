
% dlmwrite('aMat.txt', A)
H = hilb(n);
X0 = [1, 2, 3];

Ad = det(A);
if Ad == 0
   disp('Система не имеет единственного решения');
end

Kd = [];
for i =1:n
    K = A(1:i, 1:i);
    Kd = [Kd, det(K)];
    if det(K)==0
       disp('Матрица неразложима')
    end
end

sol = LDLTfun(A, X0);
solH = LDLTfun(H, X0);

c = cond(A);
ch = cond(H);
deltasol = norm(X0.'-sol)/norm(X0);
C =[];
D = [];
C = [C, c];
D = [D, deltasol];
for i =3:12
    T = 10*rand(n);
    [Q, R] = qr(T);
    V = [1,1, 10^i];
    D1 = diag(V);
    A1 = Q.'*D1*Q;
    c = cond(A1);
    C = [C, c];
    d = norm(X0.'- LDRfun(A1, X0))/norm(X0);
    D = [D, d];
end
deltasolH = norm(X0.'-solH)/norm(X0);
C = [C, ch];
D = [D, deltasolH];

figure 
loglog(C, D, '-*')
xlabel('cond')
ylabel('||dX||/||X||')