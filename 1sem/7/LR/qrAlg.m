%%Хорошо отделимые сч
Cr = dlmread('aMatEig.txt');
d = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1];
D = diag(d);
A = Cr^-1*D*Cr;
[e11,EV1] = LRmethod(A);

(A-10*eye(10))*EV1(:,1)