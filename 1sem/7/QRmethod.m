function [e, V, t2, t1] = QRmethod(A, H, q, er, eps0)
Ao = A;
Ho = H;
K = [];
n = length(A);
V= [];
Err = zeros(n,1);

    eps = 10^-10;
    k = 0;
    k1 =0;
    A = Ao;
    H = Ho;
    tic
    while max(tril(A, -1), [], 'all')>eps
        [Q,R] = givens(A);
        A = Q.'*A*Q;
        e1 = diag(A);
        k = k+1; 
    end
    t1 = toc;
    tic 
    while max(tril(H, -1), [], 'all')>eps
        [Q,R] = givens(H);
        H = Q.'*H*Q;
        eh1 = diag(H);
        k1 = k1+1;
    end
    t2 = toc;
e = e1;

for t = 1:n
    y = RevIt(Ho, e1(t)+0.1);
    v=q*y;
    V = [V, v];
end

% figure 
% hold on 
% grid on
% for i = 1:n
%     plot(K, Err(i,:))
% end
% title('Погрешность сч от номера итерации(плохо отделимые сч)')
% xlabel('Кол-во итераций')
% ylabel('Ошибка')
% figure
% semilogx(EPS, K, '-r')
% hold on
% grid on
% semilogx(EPS, K1, '-b')
% % for i=1:n
% %     plot(
% title('Зависимость кол-ва итераций от точности')
% legend('Изначальная матрица', 'Форма Хессенберга')
% xlabel('Eps')
% ylabel('Кол-во итераций')
end