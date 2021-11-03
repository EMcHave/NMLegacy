function [Q, R] = givensMod(A)
n = length(A);
Q = eye(n);
for i = 1:n-1
        t=A(i,i)/A(i+1,i);
        c = 1/(1+t^2)^0.5;
        s = t*c;
        q = eye(n);
        q(i,i) = s;
        q(i+1,i+1) = s;
        q(i,i+1)= c;
        q(i+1,i) = -c;
        A = q*A;
        Q = q*Q;
end
Q = Q.';
R = A;
end