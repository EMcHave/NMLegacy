A = dlmread('qMat.txt');
n = length(A);
X0 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; %столбец - точное решение
B = A*X0.';

H = hilb(n);

Ad = det(A);
if Ad == 0
   disp('Система не имеет единственного решения');
end


C =[]; %массив чисел обусловленности
D = []; %массив относительных погрешностей
Dgiv = [];

D2 = [];

dX = []; % массив погрешностей решения при возмущении столбца свободных членов
dXgiv = [];
dB = []; % массив погрешностей столбца свободных членов


for i =3:11
    T = 10*rand(n);
    [Q, R] = qr(T);
    
    V = [1,1,1,1,1,1,1,1,1,10^i]; %вектор с.ч
    D1 = diag(V);                 %диагональная матрица на основании вектора V
    A1 = Q.'*D1*Q;
    c = cond(A1);
    C = [C, c];
    d = norm(X0.'- gauss_column(A1, X0))/norm(X0);
    dgiv = norm(X0.' - givensQRfun(A1, X0))/norm(X0);
    D = [D, d];
    Dgiv = [Dgiv, dgiv];
    d2 = norm(X0.'- gauss(A1, X0))/norm(X0);
    D2 = [D2, d2];
end

for i = 10:-1:-2
    %In general, you can generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1).
    dist = -20 + 40.*rand(n, 1); %вектор возмущений - от -20 до 20
    Bdist = B + (10^-i)*dist;
    dB = [dB, norm(B - Bdist)/norm(B)];
    dX = [dX, norm(X0.' - gauss_columnForB(A, Bdist))/norm(X0)];
    dXgiv = [dXgiv, norm(X0.' - givensQRfunForB(A, Bdist))/norm(X0)];
end

% сравнить методы Гаусса и мод. придумать матрицу, чтобы была видна разница
% (для 1 исследования)
figure 

loglog(C, D, '-*')
hold on
loglog(C, D2, 'r')
loglog(C, Dgiv, 'g')
grid on
xlabel('cond')
ylabel('||dX||/||X||')

figure 
loglog(dB, dX, '-*')
hold on
loglog(dB, dXgiv, 'r')
grid on
xlabel('||dB||/||B||')
ylabel('||dX||/||X||')