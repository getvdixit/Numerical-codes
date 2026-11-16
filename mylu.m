% This function finds the lu-decomposition of the matrix A.

function [L,U]=mylu(A)

[m,n] = size(A);
L = eye(m);
U = A;
for i = 1:m-1
   if U(i,i) == 0 error('Pivot is 0'); return; end
   for k = i+1:m
       L(k,i) = U(k,i)/U(i,i);
	   for j = i+1:n
	       U(k,j) = U(k,j) - L(k,i)* U(i,j);
	   end
	   for j = 1:i
	       U(k,j) = 0;
	   end
   end
end
