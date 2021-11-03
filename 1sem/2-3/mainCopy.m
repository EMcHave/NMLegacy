A = dlmread('aMat.txt');
n = length(A);
% dlmwrite('aMat.txt', A)
H = hilb(n);
X0 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; %столбец - точное решение
B = A*X0.';
sol = gauss_column(A, X0);

Ad = det(A);
if Ad == 0
   disp('Система не имеет единственного решения');
end


C =[]; %массив чисел обусловленности
D = []; %массив относительных погрешностей

dX = []; % массив погрешностей решения при возмущении столбца свободных членов
dB = []; % массив погрешностей столбца свободных членов


for i =3:11
    T = 10*rand(n);
    [Q, R] = qr(T);
    %число обусловленности пропорционально отношению максиального с.ч и
    %минимального

    %чем больше i - тем больше отношение (10^i)/1 - тем больше число
    %обусловенности - тем более плохо обусловена матрица
    V = [1,1,1,1,1,1,1,1,1,10^i]; %вектор с.ч
    D1 = diag(V);                 %диагональная матрица на основании вектора V
    A1 = Q.'*D1*Q;
    c = cond(A1);
    C = [C, c];
    d = norm(X0.'- gauss_column(A1, X0))/norm(X0);
    D = [D, d];
end

for i = 10:-1:-2
    %In general, you can generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1).
    dist = -20 + 40.*rand(n, 1); %вектор возмущений - от -20 до 20
    Bdist = B + (10^-i)*dist;
    dB = [dB, norm(B - Bdist)/norm(B)];
    dX = [dX, norm(X0.' - gauss_columnForB(A, Bdist))/norm(X0)];
end


figure 

loglog(C, D, '-*')
hold on
grid on
xlabel('cond')
ylabel('||dX||/||X||')

figure 
loglog(dB, dX)
hold on
grid on
xlabel('||dB||/||B||')
ylabel('||dX||/||X||')