function [newEr] = methodsEr(x10, x20, eps0, f, g, a1, b1, a2, b2)
syms x
df = diff(f(x));
dg = diff(g(x));
ddf = diff(df);
ddg = diff(dg);

M20 = double(subs(ddf, b1));
m10 = double(subs(df, a1));

M21 = double(subs(ddg, b2));
m11 = double(subs(dg, a2));

k1 = 0;
k2 = 0;
n1 = 0;
n2 = 0;

Ke1 = [];
Ke2 = [];
Ne1 = [];
Ne2 = [];

Err1 = [];
Err2 = [];
Nerr1 = [];
Nerr2 = [];

rt1 = fzero(f, x10);
rt2 = fzero(g, x20);

x11 = x10 - f(x10)/double(subs(df,x10));
x21 = x20 - g(x20)/double(subs(dg,x20));
while M20/(2*m10)*(abs(x10 - x11))^2 > eps0
        t = x11;
        x10 = x11;
        x11 = t - (f(t)/double(subs(df,t)));
        err1 = abs(rt1 - x11);
        Err1 = [Err1, err1];
        Ke1 = [Ke1, k1];
        k1 = k1 + 1;
end

while M21/(2*m11)*(abs(x20 - x21))^2 > eps0
        t = x21;
        x20 = x21;
        x21 = t - (g(t)/double(subs(dg,t)));
        err2 = abs(rt2 - x21);
        Err2 = [Err2, err2];
        Ke2 = [Ke2, k2];
        k2 = k2 + 1;
end
while abs(a1-b1) > 2*eps0
     y = (a1+b1)/2;
     if f(a1)*f(y)<0 
         b1 = y;
     else
         a1 = y;
     end
     nerr1 = abs(rt1 - y);
     Nerr1 = [Nerr1, nerr1];
     Ne1 = [Ne1, n1];
     n1 = n1+1;
end
 
while abs(a2-b2) > 2*eps0
     y1 = (a2+b2)/2;
     if g(a2)*g(y1)<0 
         b2 = y1;
     else
         a2 = y1;
     end
     nerr2 = abs(rt2 - y1);
     Nerr2 = [Nerr2, nerr2];
     Ne2 = [Ne2, n2];
     n2 = n2+1;
 end
figure 
semilogx(Ke1, Err1, 'b')
hold on
semilogy(Ke2, Err2, 'r')
semilogy(Ne1, Nerr1, 'g')
semilogy(Ne2, Nerr2, 'm')
legend('f, newton', 'g, newton', 'f, bisection', 'g, bisection')
title('Зависимость ошибки от k')
newEr = [Ke1, Ke2, Err1, Err2];
end

