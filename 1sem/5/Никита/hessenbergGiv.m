function Hess = hessenbergGiv(A)
n = length(A);
Q = eye(n);
for i = 2:n-1
    for j = i+1:n
        c = A(i,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
        s = A(j,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
        q = eye(n);
        q(i,i) = c;
        q(j,j) = c;
        q(i,j)= -s;
        q(j,i) = s;
        A = q.'*A*q;
        Q = Q*q;
    end
end
Hess = A;
end
