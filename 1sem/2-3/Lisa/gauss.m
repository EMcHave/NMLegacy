function solution = gauss(matrix, X0)

free = matrix*X0.';
n = size(matrix, 1);
A = cat(2, matrix, free);

for j = 1:n 
    for i = (j+1):n    
        m = A(i, j) / A(j, j);
        A(i,:) = A(i,:) - m*A(j,:);
    end 
end

solution = zeros(n,1);
for s = n:-1:1
   numerator = A(s,end) ;
   for numer = 1:n
       if numer ~= s
           numerator = numerator - A(s, numer)*solution(numer ,1);
       end
   end
   solution(s,1) = numerator / A(s,s); 
end