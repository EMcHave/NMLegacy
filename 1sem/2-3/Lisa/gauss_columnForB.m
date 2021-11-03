function solution = gauss_columnForB(matrix, free)

n = size(matrix, 1);
A = cat(2, matrix, free);

for j = 1:n 
    
    %выбор ведущего элемента столбца
    [M, I] = max(abs(A(j:n,j)));
    if (I ~= j && M ~= A(j,j))
        tmp = A(j,:);
        A(j,:)  = A((j+I-1), :);
        A((j+I-1), :) = tmp;
    end
    %конеца выбора
    
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