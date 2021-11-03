f = @(x) sqrt(sin(x.^2));
a = -1;
b = 1;
epsMax = 14;
Eps = [];
K = [];
x = linspace(a, b, 500);

for i = 1:epsMax
    eps = 10^-i;
    [k, n] = runge(f, a, b, eps);
    K = [K, [eps, k, n]'];
end


I = Simpson(f,a,b,60)

figure
hold on
grid on
plot(x, f(x))
title('Function plot')
xlabel('x')
ylabel('y')

function [k, n] = runge(f,a, b, eps)
    n = 2;
    d = 1;
    k = 0;
    Ip = gauss(f, a, b, n);
    while d > eps
        In = Simpson(f,a,b,2*n);
        d = abs(Ip-In)/15;
        Ip = In;
        n = n*2;
        k = k+1;
    end
end

