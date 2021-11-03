function [RtVector] = givensQRfunForB(A, B)
n = length(A);

[Q, R] = givens(A);

% dif = eye(n)-Q*Q.'

X = [];
k = n-1;
Y = Q.'*B;
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

[RtVector] = X.';
end