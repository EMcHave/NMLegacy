A = dlmread('aMatEig.txt');
n = length(A);
xp = ones(n,1);
L1 = ones(n,1);
L2 = xp+L1;

d = [1,2,3,4,5,6,7, 8, 9, 9+10^-10];
D = diag(d);
[Q,R] = qr(A);
A1 = A^-1.'*D*A;
[V, Lam] = eig(A);
[V1, Lam1] = eig(A1);
[lm, ev]  = powerrrr(A, xp, L1, L2);  %cч и св для матрицы А (хорошо отдлимое сч)
[lm1, ev1] = powerrrr(A1, xp, L1, L2);  %cч и св для матрицы А1 (плохо отдлимое сч)

%проверка (A-lambda*E)=0
disp((A-lm*eye(n))*ev)          %для матрицы A
disp((A1-lm1*eye(n))*ev1)       %для матрицы А1