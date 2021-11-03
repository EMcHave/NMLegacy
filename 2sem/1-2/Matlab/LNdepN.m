function [Er] = LNdepN(f,a,b,ns, nf)
Er = [];
Er1 = [];
N = [];
N1 = [];
x=a:0.01:b;
for n = ns:nf 
%%  Нормальная сетка   
    h = (b-a)/n;
    Xi = []; 
    Yi = [];
    Xi = [Xi, a];
    for i=1:n
        Xi = [Xi, a+i*h];
        Yi = [Yi, f(Xi(i))];
    end
    Yi = [Yi, f(b)];
    sum=0;
    for k=1:length(Xi)
        p=1;
        for j=1:length(Xi)
            if j~=k
                c = poly(Xi(j))/(Xi(k)-Xi(j));
                p = conv(p,c);
            end
        end
        term = p*Yi(k);
        sum= sum + term;
    end
    LagPoly1 = sum;
    
    E = @(x) polyval(LagPoly1,x)-f(x);
    Err = E(x);
    Er = [Er, max(Err)];
    N = [N,n];
%%  Чебышевская сетка
    Xi1 = []; 
    Yi1 = [];
    Xi1 = [Xi1, (b+a)/2+((b-a)/2)*cos(pi/(2*n+2))];
    for i=1:n
        xi1 = (b+a)/2+((b-a)/2)*cos(((2*i+1)/(2*n+2))*pi);
        Xi1 = [Xi1, xi1];
        Yi1 = [Yi1, f(Xi1(i))];
    end
    Yi1 = [Yi1, f(b)];
    sum1=0;
    for k=1:length(Xi1)
        p1=1;
        for j=1:length(Xi)
            if j~=k
                c1 = poly(Xi1(j))/(Xi1(k)-Xi1(j));
                p1 = conv(p1,c1);
            end
        end
        term1 = p1*Yi1(k);
        sum1= sum1 + term1;
    end
    LagPoly11 = sum1;
    E1 = @(x) polyval(LagPoly11,x)-f(x);
    Err1 = E1(x);
    Er1 = [Er1, max(Err1)];
    N1 = [N1,n];
end
figure
hold on
grid on
plot(N, Er, '*-r')
plot(N1, Er1, '*-b')
title('Зависимость макс. ошибки от числа узлов(5-45 узлов)')
xlabel('N')
ylabel('Err')
legend('Р/м сетка','Чебышевская')
end