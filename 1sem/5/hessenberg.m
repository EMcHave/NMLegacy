function Hess = hessenberg(A)
n = length(A);
H = A;
sum = 0;
for m = 1:n-2
    for i=m+1:n
        sum = sum+(H(i,m))^2;
    end
    s = sign(-A(m+1,m))*sqrt(sum);
    w = H(:,m);
    for j=1:m
        w(j)=0;
    end
    w(m+1) = w(m+1)-s;
    w = w/(sqrt(2*s*(s-H(m+1,m))));
    Q = eye(n)-2*w*w.';
    H = Q*H*Q;
    sum = 0;
end
Hess = H;
end

