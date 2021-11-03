function [RtVector] = GSnmQRfun(A, X0)
X0 = X0.';
B = A*X0;
n = length(A);
Q = zeros(n);
R = zeros(n);
A1 = A;

for i = 1:n
    Q(:,i) = A(:,i);
    for j = 1:i-1
        R(j, i) = dot(Q(:,j),A(:,i));
        Q(:,i) = Q(:,i) - R(j,i)*Q(:,j);
    end
    R(i,i) = norm(Q(:,i));
    if R(i,i) == 0
        break
    end
    Q(:,i) = Q(:,i)/R(i,i);
end


X = [];
k = n-1;
while k > 0
    Y = Q.'*B;
    X(n) = Y(n)/R(n,n);
    sum1 = 0;
    for j = k+1:n
        sum1 = sum1 + R(k,j)*X(j);
    end
    X(k) = (Y(k)-sum1)/R(k,k);
    k = k-1;
    sum1 = 0;
end



[RtVector] = X.';
end