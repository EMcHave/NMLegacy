function [newErRt] = newtonErRt(x0, root, eps0, f)

syms x
df = diff(f(x));
x10 = x0;
k1 = 0;
k2 = 0;

Kr1 = [];
Kr2 = [];
Er1 = [];
Er2 = [];
C= [];
C1= [];
x1 = x0 - 2*(f(x0)/double(subs(df,x0)));
x11 = x10 - (f(x10)/double(subs(df,x10)));
while (abs(x0 - x1)) > eps0
        t = x1;
        x0 = x1;
        x1 = t - 2*(f(t)/double(subs(df,t)));
        er1 = abs(root - x1);
        Er1 = [Er1, er1];
        Kr1 = [Kr1, k1];
        k1 = k1 + 1;
end

while (abs(x10 - x11)) > eps0
        t = x11;
        x10 = x11;
        x11 = t - (f(t)/double(subs(df,t)));
        er2 = abs(root - x11);
        Er2 = [Er2, er2];
        Kr2 = [Kr2, k2];
        k2 = k2 + 1;
end
for i =1:length(Er1)-1
    C = [C, Er1(i+1)/(Er1(i))^2];
end

for i =1:length(Er2)-1
    C1 = [C1, Er2(i+1)/(Er2(i))];
end
figure 

semilogy(Kr1, Er1, 'b')
hold on
semilogy(Kr2, Er2, 'r')
legend('исп. формула для кр. корня', 'стандартная формула')
title('Зависимость ошибки от k (случай кратного корня)')
newErRt = [C, C1];
end
