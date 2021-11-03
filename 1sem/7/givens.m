function [Q, R] = givens(A)
n = length(A);
Q = eye(n);
for i = 1:n
    for j = i+1:n
        c = A(i,i)/((A(i,i))^2 + (A(j,i)^2))^0.5;
        s = A(j,i)/((A(i,i))^2 + (A(j,i)^2))^0.5;
        q = eye(n);
        q(i,i) = c;
        q(j,j) = c;
        q(j,i)= -s;
        q(i,j) = s;
        A = q*A;
        Q = q*Q;
    end
end
R = A;
Q = Q.';
end