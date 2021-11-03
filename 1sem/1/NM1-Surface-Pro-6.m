X=linspace(-100,100,1000);
eps = 1e-3;
g = @(X) 5.^X - 6.*X - 7;
f1 =@(X) (X-sqrt(2)).^2.*(X.^3-7.*X-13);

syms x
df = diff(f1(x));
dg = diff(g(x));
ddf = diff(df);
ddg = diff(dg);


solNew = (methods(3, 3.1, 2, 1.9  eps, f1, g, 3, 4, 1.5, 3));
disp(methodsEr(4, 3, 1e-15, f1, g, 3, 4, 1.5, 2.5));
disp(newtonErRt(0.4,solNew(:, 1), 1e-15, f1));



figure

subplot(2,1,1)
hold on 
grid on
axis([-4, 5, -60, 30])
plot(X, f1(X))
plot(solNew(:,1), 0, 'ok')
plot(solNew(:,2), 0, 'ok')
title('График f1')


subplot(2,1,2)
hold on
grid on
axis([-5, 5, -10, 10])
plot(X, g(X))
plot(solNew(:,3), 0, 'ok')
title('График g')

figure
subplot(2,1,1)
hold on
grid on
axis([3, 4, -50, 500])
plot(X, double(subs(df, X)), 'r')
plot(X, double(subs(ddf, X)), 'b')
title('Df, ddf')

subplot(2,1,2)
hold on
grid on
axis([1.5, 2.5, 0, 300])
plot(X, double(subs(dg, X)), 'r')
plot(X, double(subs(ddg, X)), 'b')
title('Dg, ddg')


