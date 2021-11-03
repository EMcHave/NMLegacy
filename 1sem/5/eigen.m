A = dlmread('aMatEig.txt');
n = length(A);
A = A(1:6,1:6);
H = hessenberg(A);
x = -10:0.05:35;
%%  Нормально отделимые С.Ч
Char = DetH(H);
p = sym2poly(Char);
e = roots(p);
e = e(imag(e)==0);
syms t
ChP = @(x) double(subs(Char, t, x));
er = eig(A);
%% Плохо отделимые С.Ч
[Q, R] = qr(A);
d = [10, 5, 5+10^-7, 3, 3+10^-6, 1];
D = diag(d);
A1 = Q.'*D*Q;
H1 = hessenberg(A1);
Char1 = DetH(H1);
p1 = sym2poly(Char1);
e1 = roots(p1);
er1 = d;
%% Большая матрица
% [Q, R] = qr(rand(25));
% db = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25];
% Db = diag(db);
% A1b = Q.'*Db*Q;
% H1b = hessenberg(A1b);
% Char1b = DetH(H1b);
% p1b = sym2poly(Char1b);
% eb1 = roots(p1b);
% erb1 = db;
%%
figure
hold on
grid on
axis([-10, 35, -120000, 50000])
plot(x, double(subs(Char, x)))
plot(e, 0, '*')

figure
hold on
grid on
axis([-1, 15, -0.5, 0.5])
plot(x, double(subs(Char1, x)))
plot(e1, 0, '*')
