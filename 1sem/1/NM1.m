X=linspace(-100,100,1000);
eps = 1e-3;
g = @(X) 5.^X-2-exp(-2.*X);
f1 =@(X) (X-sqrt(2)).^2.*(X.^3-7.*X-13);

syms x
df = diff(f1(x));
dg = diff(g(x));
ddf = diff(df);
ddg = diff(dg);


sol = (methods(4, 3.9, 1, 0.9, eps, f1, g, 3, 4, 0.4, 1.5));
disp(methodsEr(4, 3.9, 1, 0.9, 10e-15, f1, g, 3, 4, 0.4, 1.5));



figure

subplot(2,1,1)
hold on 
grid on
axis([-4, 5, -60, 30])
plot(X, f1(X))
plot(sol(:,1), 0, 'ok')
title('График f1')


subplot(2,1,2)
hold on
grid on
axis([-5, 5, -10, 10])
plot(X, g(X))
plot(sol(:,2), 0, 'ok')
title('График g')

figure
subplot(2,1,1)
hold on
grid on
axis([3, 4, -50, 500])
plot(X, double(subs(df, X)), 'b')
plot(X, double(subs(ddf, X)), 'r')
legend('df', 'ddf')
title('Df, ddf')

subplot(2,1,2)
hold on
grid on
axis([0.4, 1.5, 0, 50])
plot(X, double(subs(dg, X)), 'b')
plot(X, double(subs(ddg, X)), 'r')
legend('dg', 'ddg')
title('Dg, ddg')


