m = 6;
p = floor((m+1)/2);

A = zeros(2*p+1, 2*p+1);
for i=1:(2*p+1)
   for j=1:(2*p+1)
       A(j,i) = (i-p-1)^(j-1);
   end
end
A;

b = zeros([2*p+1, 1]);
b(m+1) = factorial(m);
b;
A\b