

F = [1,1;1, 1.0001]
fc = cond(F)

A = hilb(8)

X = ones(8, 1)

B = A*X

X1 = A\B

dx = X-X1;


dx = norm(dx)/norm(X)
cond(A)