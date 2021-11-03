function [new] = methods(x10, x11, x20, x21, eps0, f, g, a1, b1, a2, b2)


syms x
df = diff(f(x));
dg = diff(g(x));
ddf = diff(df);
ddg = diff(dg);


i = 0;
k = 0;
k1 = 0;
k2 = 0;
n = 0;
n1 = 0;

K1 = [];
K2 = [];
N1 = [];
N2 = [];
eps = eps0;
EPS = [];


M20 = double(subs(ddf, b1));
m10 = double(subs(df, a1));

M21 = double(subs(ddg, b2));
m11 = double(subs(dg, a2));

x12 = x11 - (f(x11)*(x11-x10))/(f(x11)-f(x10));
x22 = x21 - g(x21)*(x21-x20)/(g(x21)-g(x20));
while eps > 10e-15

    while M20/(2*m10)*(abs(x10 - x11)) > eps
        t = x10;
        h = x11;
        x10 = x11;
        x12 = h - f(h)*(h-t)/(f(h)-f(t));
        x11 = x12;
        k1 = k1 + 1;
    end

    while M21/(2*m11)*(abs(x20 - x21)) > eps
        t = x20;
        h = x21;
        x20 = x21;
        x22 = h - g(h)*(h-t)/(g(h)-g(t));
        x21 = x22;
        k2 = k2 + 1;
    end
    
    
    
    while abs(a1-b1) > 2*eps
     y = (a1+b1)/2;
        if f(b1)*f(y)<0 
            b1 = y;
        else
            a1 = y;
        end
    n = n+1;
    end
 
    while abs(a2-b2) > 2*eps
        y1 = (a2+b2)/2;
        if g(a2)*g(y1)<0 
            b2 = y1;
        else
            a2 = y1;
        end
    n1 = n1+1;
    end
    
    N1 = [N1, n];
    N2 = [N2, n1];
    K1 = [K1, k1];
    K2 = [K2, k2];
    EPS = [EPS, eps];
    eps = 10^(-3-i);
    i = i+1;
end
new = [x12, x22];

figure 
semilogx(EPS, K1, 'b')
hold on
semilogx(EPS, K2, 'r')
semilogx(EPS, N1, 'g')
semilogx(EPS, N2, 'y')
legend('f, newton', 'g, newton', 'f, bisection', 'g, bisection')
title ('Зависимость k от eps')

end

