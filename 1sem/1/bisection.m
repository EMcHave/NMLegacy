function [bis] = bisection(a, b, a1, b1, f, g, eps)
n = 0;
n1 = 0;
 while abs(a-b) > 2*eps
     c = (a+b)/2;
     if f(a)*f(c)<0 
         b = c;
     else
         a = c;
     end
     n = n+1;
 end
 
  while abs(a1-b1) > 2*eps
     c1 = (a1+b1)/2;
     if g(a1)*g(c1)<0 
         b1 = c1;
     else
         a1 = c1;
     end
     n1 = n1+1;
  end
  bis = [c, c1];
end

