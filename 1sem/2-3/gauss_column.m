function solution = gauss_column(matrix, X0)

%A*X0 = B
%matrix*X0 = free

free = matrix*X0.'; %столбец свободных членов на основании уже известного решения
n = size(matrix, 1);
A = cat(2, matrix, free);

for j = 1:n 
    
    %choosing of leading element of column
    [M, I] = max(abs(A(j:n,j)));
    if (I ~= j && M ~= A(j,j))
        tmp = A(j,:);
        A(j,:)  = A((j+I-1), :);
        A((j+I-1), :) = tmp;
    end
    % end of the choosing
    
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