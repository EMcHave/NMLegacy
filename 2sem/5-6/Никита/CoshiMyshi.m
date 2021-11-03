d2y = @(x, y, dy) dy*sin(x) - y*cos(x) - cos(x) + 1;
y_exact = @(x) cos(x);
y0 = 1;
dy0 = 0;
a = 0;
b = pi/2;
X = linspace(a, b, 200);

[x, y, z] = coshimyshi(d2y, y0, dy0, a, b, 1);
[Eps, EpsR, N] = runge(y_exact, 10, d2y, y0, dy0, a, b);


figure
semilogx(Eps, N)
grid on
xlabel('Eps')
ylabel('k')
title('Количество итераций от точности')

figure
loglog(Eps, EpsR)
grid on
xlabel('Eps')
ylabel('Абс Погр')
title('Абсолютная погрешность от заданной точности')

function [Eps, EPSREAL, N] = runge(y_exact, epsMax, d2y, y0, dy0, a, b)
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
            for j =2:length(yp)
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