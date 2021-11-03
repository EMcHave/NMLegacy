function [RtVector] = LUfun(A, X0)
X0 = X0.';
B = A*X0;
n = length(A);
L = eye(n);
D = eye(n);
for m = 1:n
    S1 = 0;
    S2 = 0;
    for j = m:n 
        for g = 1:m-1
            S1 = S1+L(m, g)*R(g, j);
        end
        R(m, j) = A(m, j)-S1;
        S1 = 0;
    end
    
    for i = m+1:n
        for t = 1:m-1
            S2 = S2+L(i, t)*R(t, m);
        end
        L(i, m) = (A(i, m)-S2)/R(m,m);
        S2=0;
    end
end 

X = [];
Y = [];
for i = 2:n
    sum = 0;
    Y(1) = B(1);
    for j = 1:i-1
        sum = sum + L(i, j)*Y(j);
    end
    Y(i) = B(i) - sum;
    sum = 0;
end

k = n-1;
while k > 0
    X(n) = Y(n)/R(n,n);
    sum1 = 0;
    for j = k+1:n
        sum1 = sum1 + R(k,j)*X(j);
    end
    X(k) = (Y(k)-sum1)/R(k,k);
    k = k-1;
    sum1 = 0;
end
X.'

RtVector = X.';
end

