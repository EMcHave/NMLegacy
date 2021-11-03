function In = Simpson(f, a, b, n)
    In = 0;
    h = (b-a)/n;
    X = [];
    for i = 0:n
        X = [X, a+i*h];
    end
    for i = 1:(n/2)
        In = In + h/3*( f(X(2*i-1)) + 4*f( (X(2*i) )) + f(X(2*i+1)) );
    end
end
