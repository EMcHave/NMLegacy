function [LagPoly1, Xi, Yi] = lagrangeCheb(f,a,b,n)

Xi = []; 
Yi = [];
Xi = [Xi, (b+a)/2+((b-a)/2)*cos(pi/(2*n+2))];
for i=1:n
    xi = (b+a)/2+((b-a)/2)*cos(((2*i+1)/(2*n+2))*pi);
    Xi = [Xi, xi];
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
end
