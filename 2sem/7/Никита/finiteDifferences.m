yR = @(x) cos(x);
a = 0;
b = pi/2;
ya = 1; 
yb = yR(b);
p = @(x) -sin(x);
q = @(x) cos(x);
f = @(x) 1 - cos(x);
[x, y, n] = runge(3, p, q, f, ya, yb, a, b);
X = linspace(a, b, 100);

[x1, y1] = finDif(a, b, 16, p, q, f, ya, yb); %sol10
[x2, y2] = finDif(a, b, 8, p, q, f, ya, yb); %sol5
Eps = [];
Epsr = [];
table = [];
for i = 1:length(y2)
    E = abs(y2(i) - y1(2*i-1))/3;
    Er = abs(yR(x1(2*i-1)) - y1(2*i-1));
    v = [x1(2*i-1); E; Er];
    table = [table, v];
    Eps = [Eps, E];
    Epsr = [Epsr, Er];
end

tableDist = [];
[xud, yud] = finDif(a, b, 6, p, q, f, ya, yb);
for i = 1:6
   r = 0.1*rand;
   yAdist = ya + r;
   [xd, yd] = finDif(a, b, 6, p, q, f, yAdist, yb);
   diff = [];
   diff = [diff; r];
   for j = 1:length(yd)
      diff = [diff; abs(yud(j) - yd(j))]; 
   end
   tableDist = [tableDist, diff];
end    




figure
hold on
grid on
plot(x, Eps)
plot(x, Epsr)
xlabel('x')
ylabel('Погрешность')
legend('По правилу Рунге','Действительная')


function [xn, yn, n] = runge(eps, p, q, f, ya, yb, a, b)
        eps = 10^-eps;
        n = 2;
        d = 1;
        [xp, yp] = finDif(a, b, n, p, q, f, ya, yb);
        while d > eps
            D = [];
            [xn, yn] = finDif(a, b, 2*n, p, q, f, ya, yb);
            for j =1:length(yp)-1
                d = abs(yn(2*j-1) - yp(j))/3;
                D = [D, d];
            end
            yp = yn;
            n = 2*n;
            d = max(D);
        end
        
end

function [x, y] = finDif(a, b, n, p, q, f, ya, yb) 
    h = (b-a)/n;
    x = [];
    y = [];
    %delta, lambda - прогоночные коэффициенты    
    lambda = [];
    delta = [];
    for i = 0:n
        x = [x, a+i*h];
        y = [y, 0];
        lamdda = [lambda, 0];
        delta = [delta, 0];
    end
    y(1) = ya;
    y(n+1) = yb;
    lambda(1) = ya;
    % Прямой ход метода прогонки     
    for i = 2:n
        B = 1 - h/2 * p(x(i));
        C = h*h*q(x(i)) - 2;
        D = 1 + h/2 * p(x(i));
        R = h*h*f(x(i));
        delta(i) = -D/(B*delta(i-1)+C);
        lambda(i) = (R - B*lambda(i-1))/(B*delta(i-1)+C);
    end
    % обратный ход     
    for j = n:-1:2
        y(j) = delta(j) * y(j+1) + lambda(j);
    end
end