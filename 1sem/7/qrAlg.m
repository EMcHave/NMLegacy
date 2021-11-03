%%Хорошо отделимые сч
Cr = dlmread('aMatEig.txt');
d = [10,9,8,7,6,5,4,3,2,1];
D = diag(d);
A = Cr^-1*D*Cr;
eps0 = 10^-12;
[H, q] = hessenbergGiv(A);
[e1,EV1] = QRmethod(A,H,q, d.', eps);
%%Плохо отделимые сч
d1 = [9+10^-11, 9, 8 , 6+10^-11, 6 , 5, 4, 2+10^-10, 2,1];
D1 = diag(d1);
A1 = Cr^-1*D1*Cr;
[H1,q1] = hessenberg(A1);
[e2, EV2, tA, tH] = QRmethod(A1,H1,q1, d1.');

dV=[];
dV1 = [];
for i = 1:10
  dV = [dV, (A-e1(i,i)*eye(10))*EV1(:,i)];
  dV1 = [dV1, (A-e2(i,i)*eye(10))*EV2(:,i)];
end

syms y x a
p = x^3/3+y*x^2-5*x+y^2/3-5*y - 10
factor(p)