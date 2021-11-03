function [RtVector] = LDLTfun(A, X0)
X0 = X0.';
B = A*X0;
n = length(A);
L = eye(n);
D = eye(n);
for m = 1:n
    S = 0;
    S2 = 0;
    R = L.';
    for k = 1:m-1
        S = S + L(m, k)*D(k, k)*R(k,m);
    end
    D(m,m) = A(m,m) - S;
    
    for i = m+1:n
        for k = 1:m-1
            S2 = S2+L(i, k)*D(k, k)*R(k, m);
        end
        L(i, m) = (A(i, m)-S2)/D(m,m);
        S2 = 0;
    end
end
L*D*L.'
X = [];
Y = [];
Z = [];
for i = 2:n
    sum = 0;
    Z(1) = B(1);
    for j = 1:i-1
        sum = sum + L(i, j)*Z(j);
    end
    Z(i) = B(i) - sum;
    sum = 0;
end
Z(n)
for i = 1:n
    Y(i) = Z(i)/D(i,i);
end

k = n-1;
while k > 0
    X(n) = Y(n);
    sum1 = 0;
    for j = k+1:n
        sum1 = sum1 + R(k,j)*X(j);
    end
    X(k) = (Y(k)-sum1);
    k = k-1;
    sum1 = 0;
end

RtVector = X.';
end

