A = dlmread('aMat.txt');
n = length(A);
% dlmwrite('aMat.txt', A)
H = hilb(n);
X0 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

Ad = det(A);
if Ad == 0
   disp('Система не имеет единственного решения');
end

% Kd = [];
% for i =1:n
%     K = A(1:i, 1:i);
%     Kd = [Kd, det(K)];
%     if det(K)==0
%        disp('Матрица неразложима')
%     end
% end
% 
% sol = LDRfun(A, X0);
solQ = GSQRfun(A, X0);
hilnmQ = GSnmQRfun(H, X0);
hilQ = GSQRfun(H, X0);
c = cond(A);
ch = cond(H);
% deltasol = norm(X0.'-sol)/norm(X0);

C =[];
% D = [];
DQ = [];
DQnm = [];
C = [C, c];
% D = [D, deltasol];
DQ = [DQ, norm(X0.'-solQ)/norm(X0)];
DQnm = [DQnm, norm(X0.'-GSnmQRfun(A, X0))/norm(X0)];
for i =3:11
    T = 10*rand(n);
    [Q, R] = qr(T);
    V = [1,1,1,1,1,1,1,1,1,10^i];
    D1 = diag(V);
    A1 = Q.'*D1*Q;
    c = cond(A1);
    C = [C, c];
%     d = norm(X0.'- LDRfun(A1, X0))/norm(X0);
    d1 = norm(X0.'- GSQRfun(A1, X0))/norm(X0);
%     D = [D, d];
    DQ = [DQ, d1];
    DQnm = [DQnm, norm(X0.'-GSnmQRfun(A1, X0))/norm(X0)];
end
% deltasolH = norm(X0.'-LDRfun(H, X0))/norm(X0);
% deltasolHQ = norm(X0.'-GSQRfun(H, X0))/norm(X0);

% C = [C, ch];
% D = [D, deltasolH];
% DQ = [DQ, deltasolHQ];
% DQnm = [DQnm, norm(X0.'-GSnmQRfun(H, X0))/norm(X0)];

figure 

% loglog(C, D, '-*')

loglog(C, DQ, '-o')         %модифицирванный ГШ
hold on
loglog(C, DQnm, '-k')       %немодифицированный
xlabel('cond')
ylabel('||dX||/||X||')