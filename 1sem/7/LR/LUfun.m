function [L, R] = LUfun(A)
n = length(A);
L = eye(n);
R = eye(n);
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
end

