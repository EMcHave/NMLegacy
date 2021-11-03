d2y = @(x, y, dy) dy*sin(x) - y*cos(x) - cos(x) + 1;
y_exact = @(x) cos(x);
y0 = 1;
dy0 = 0;
a = 0;
b = pi/2;
X = linspace(a, b, 200);

[x, y, z] = adams(d2y, y0, dy0, a, b, 2);
[Eps, EpsR, N] = rungeAdams(y_exact, 10, d2y, y0, dy0, a, b);
[EpsC, EpsRC, NC] = rungeEuler(y_exact, 10, d2y, y0, dy0, a, b);


figure
semilogx(Eps, N)
hold on
semilogx(EpsC, NC)
grid on
xlabel('Eps')
ylabel('k')
title('Количество итераций от точности')
legend('Адамс','Эйлер-Коши')

figure
loglog(Eps, EpsR)
hold on
loglog(EpsC, EpsRC)
grid on
xlabel('Eps')
ylabel('Абс Погр')
title('Абсолютная погрешность от заданной точности')
legend('Адамс','Эйлер-Коши')

% figure
% grid on
% hold on
% plot(X, y_exact(X))
% plot(x, y, 'b')


function [Eps, EPSREAL, N] = rungeAdams(y_exact, epsMax, d2y, y0, dy0, a, b)
    eps = 0.1;
    EPSREAL = [];
    Eps = [];
    N = [];
    for i = 1:epsMax
        eps = 10^-i;
        n = 2;
        d = 1;
        k=0;
        [xp, yp, zp] = adams(d2y, y0, dy0, a, b, n);
        while d > eps
            D = [];
            [xn, yn, zn] = adams(d2y, y0, dy0, a, b, 2*n);
            for j =1:length(yp)-1
                d = abs(yn(2*j-1) - yp(j))/3;
                D = [D, d];
            end
            yp = yn;
            n = 2*n;
            d = max(D);
            k=k+1;
        end
        
        N = [N, k];
        Eps = [Eps, eps];
        EPSREAL = [EPSREAL, RealEps(y_exact, xn, yn)];
    end
end

function [Eps, EPSREAL, N] = rungeEuler(y_exact, epsMax, d2y, y0, dy0, a, b)
    eps = 0.1;
    EPSREAL = [];
    Eps = [];
    N = [];
    for i = 1:epsMax
        eps = 10^-i;
        n = 2;
        d = 1;
        k=0;
        [xp, yp, zp] = coshimyshi(d2y, y0, dy0, a, b, n);
        while d > eps
            D = [];
            [xn, yn, zn] = coshimyshi(d2y, y0, dy0, a, b, 2*n);
            for j =1:length(yp)-1
                d = abs(yn(2*j-1) - yp(j))/3;
                D = [D, d];
            end
            yp = yn;
            n = 2*n;
            d = max(D);
            k=k+1;
        end
        
        N = [N, k];
        Eps = [Eps, eps];
        EPSREAL = [EPSREAL, RealEps(y_exact, xn, yn)];
    end
end

function trueeps = RealEps(y_exact, x, y)
    ep = [];
    for i=1:length(y)
        ep = [ep, abs(y_exact(x(i))-y(i))];
    end
    trueeps = max(ep);        
end


function [x,y, z] = coshimyshiStart(d2y, y0, dy0, a, b, n)
    x = [];
    x = [x,a];
    y = [];
    y = [y, y0];
    z = [];
    z = [z, dy0];
    
    h = (b-a)/n;
    for i = 1:n
        x = [x, a+i*h];
    end
    
    for i = 1:2
        y_rough = y(i) + h*z(i);
        z_rough = z(i) + h*d2y(x(i), y(i), z(i));
        y = [y, y(i) + h/2*(z(i) + z_rough)];
        z = [z,z(i) + h/2*(d2y(x(i),y(i), z(i)) + d2y(x(i+1), y_rough, z_rough))];
    end
end

function [x,y, z] = coshimyshi(d2y, y0, dy0, a, b, n)
    x = [];
    x = [x,a];
    y = [];
    y = [y, y0];
    z = [];
    z = [z, dy0];
    
    h = (b-a)/n;
    for i = 1:n
        x = [x, a+i*h];
    end
    
    for i = 1:n
        y_rough = y(i) + h*z(i);
        z_rough = z(i) + h*d2y(x(i), y(i), z(i));
        y = [y, y(i) + h/2*(z(i) + z_rough)];
        z = [z,z(i) + h/2*(d2y(x(i),y(i), z(i)) + d2y(x(i+1), y_rough, z_rough))];
    end
end

function [x, y, z] = adams(d2y, y0, dy0, a, b, n)
    h = (b-a)/n;
    x = [];
    x = [x, a];
    y = [];
    
    z = [];
    for i = 1:n
        x = [x, a+i*h];
    end
    [xC, yC, zC] = coshimyshiStart(d2y, y0, dy0, a, b, n);
    for i =1:2
        y = [y, yC(i)];
        z = [z, zC(i)];
    end
    for i =2:n
        yB = y(i) + h/2*(3*z(i) - z(i-1));
        zB = z(i) + h/2*(3*d2y(x(i),y(i),z(i)) - d2y(x(i-1), y(i-1), z(i-1)));
        y = [y, y(i) + h/2*(zB+z(i))];
        z = [z, z(i) + h/2*(d2y(x(i+1), yB, zB) + d2y(x(i), y(i), z(i)))];
    end
end
