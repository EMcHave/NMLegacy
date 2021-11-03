function [RtVector] = GSQRfun(A, X0)

X0 = X0.';
B = A*X0;
n = length(A);
Q = zeros(n);
R = zeros(n);
A1 = A;

for i = 1:n
    Q(:,i) = A1(:,i);
    for j = 1:i-1
        R(j, i) = dot(Q(:,j),A1(:,i));
        A1(:,i) = A1(:,i) - R(j,i)*Q(:,j);
    end
    R(i,i) = norm(A1(:,i));
    if R(i,i) == 0
        break
    end
    Q(:,i) = A1(:,i)/R(i,i);
end

% for i = 1:n
%     R(i,i) = norm(A1(:,i));
%     Q(:,i) = A1(:,i)/(R(i,i));
%     for j = i+1:n
%         R(i,j) = dot(Q(:,i),A1(:,j));
%         A1(:,j) = A1(:,j) - R(i,j)*Q(:,i);
%     end
% end

% [Q, R] = qr(A);

eye(n)-Q*Q.'

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