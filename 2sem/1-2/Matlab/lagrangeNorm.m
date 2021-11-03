function [LagPoly1, Xi, Yi] = lagrangeNorm(f,a,b,n)
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
end

