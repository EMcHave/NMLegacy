function ev = RevIt(A, l)
n = length(A);
ev1 = ones(n,1);
eps = 10^-12;
err = 10;
np = norm(ev1,inf);
while err>eps
    ev2 = LDRfun((A-l*eye(n)), ev1);
    [row, col] = find(abs(ev2) == norm(ev2,inf));
    nn = ev2(row, col);
    ev1 = ev2/nn;
    err = abs(np-nn);
    np = nn;
end
ev = ev1;
end



