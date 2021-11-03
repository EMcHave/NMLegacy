f = @(x) sqrt(sin(x.^2));

a = -1;
b = sqrt(pi);
[pn, Xi, Yi] = lagrangeNorm(f,a,b,25);
pol = pn;
[pchel, XiCh, YiCh] = lagrangeCheb(f,a,b,15);
Er = LNdepN(f,a,b,5, 45);
x = linspace(a, b, 500);
pn = polyval(pn, x);
pchel = polyval(pchel,x);


figure
hold on
grid on
plot(x,pn, 'r')
plot(x,pchel, 'k')
plot(x, f(x), 'b')
plot(a,0, '*')
plot(b,0, '*')
title('Interpolated function')
xlabel('x')
ylabel('y')
legend('Р/м сетка','Чебышевская','f(x)')


