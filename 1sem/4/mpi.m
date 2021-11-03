A = dlmread('aMat.txt');
n = length(A);
H = hilb(10);
% dlmwrite('aMat.txt', A)
X0 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
eps = 10^-15;

 
Keps = JacobiEps(A, X0);
Kdet = JacobiDet(A, X0, eps);
sol = Jacobi(A, X0, eps);
