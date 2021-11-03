function [newEr] = methodsEr(x10, x11, x20, x21, eps0, f, g, a1, b1, a2, b2)
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
C = [];
C1 = [];
rt1 = fzero(f, x10);
rt2 = fzero(g, x20);

x12 = x11 - (f(x11)*(x11-x10))/(f(x11)-f(x10));
x22 = x21 - g(x21)*(x21-x20)/(g(x21)-g(x20));
while M20/(2*m10)*(abs(x10 - x11))> eps0
        t = x10;
        h = x11;
        x10 = x11;
        x12 = h - f(h)*(h-t)/(f(h)-f(t));
        x11 = x12;
        err1 = abs(rt1 - x12);
        Err1 = [Err1, err1];
        Ke1 = [Ke1, k1];
        k1 = k1 + 1;
end

while M21/(2*m11)*(abs(x20 - x21))> eps0
        t = x20;
        h = x21;
        x20 = x21;
        x22 = h - g(h)*(h-t)/(g(h)-g(t));
        x21 = x22;
        err2 = abs(rt2 - x22);
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
for i = 1:(length(Err1)-1)
    C = [C, Err1(i+1)/(Err1(i))^2];
end
   
for i = 1:(length(Nerr1)-1)
    C1 = [C1, Nerr1(i+1)/(Nerr1(i))];
end

figure 
semilogy(Ke1, Err1, 'b')
hold on
semilogy(Ke2, Err2, 'r')
semilogy(Ne1, Nerr1, 'g')
semilogy(Ne2, Nerr2, 'm')
legend('f, newton', 'g, newton', 'f, bisection', 'g, bisection')
title('Зависимость ошибки от k')
newEr = Ke1;
end

