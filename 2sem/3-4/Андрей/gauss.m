function In = gauss(f,a,b, n)
    In = 0;
    h = (b-a)/n;
    X = [];
    for i = 0:n
        X = [X, a+i*h];
    end

    for i = 1:n
        In = In + gaussAlg(f,X(i), X(i+1));
    end
end

function integral = gaussAlg(f,a,b)
    integral = 0;
    xi = [-1/sqrt(3), 1/sqrt(3)];
    weighs = [1,1];
    for i = 1:2
        integral = integral + weighs(i)*f((a+b)/2 + (b-a)/2*xi(i));
    end
    integral = integral*(b-a)/2;
end
