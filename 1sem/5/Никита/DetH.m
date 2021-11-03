function Det = DetH(H)
syms t
D = [];
D = [D,1];
D = [D,(H(1,1)-t)];
n = length(H);
m = 0;
for i = 2:n
    d = 0;
    B = H(1:i, 1:i);
    K = diag(B, -1);
    for j = 1:i-1
        k=1;
        if j == 1
           k = K(end);
        else
            for v = 1:j
                k=k*K(i-v);
            end
        end
        d = d+((-1)^j)*B(i-j, i)*k*D(i-j);
    end 
    d = (B(i,i)-t)*D(i) + d;
    D = [D, d];
end
d = simplify(d);
Det = d;
end
