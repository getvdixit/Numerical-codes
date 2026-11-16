% This function calculates solution by gauss method.
%Here A is m by n(m=n) matrix and b is n by 1 matrix.


function x=gauss(A,b)
[L,U]=mylu(A);
[m,n]=size(A);

%forward substitution  L*z=b
z(1)=b(1)/L(1,1);
for i=2:n
   alpha=0;
   for j=1:(i-1)
      alpha=alpha+z(j)*L(i,j);
   end
   z(i)=(b(i)-alpha)/L(i,i);
end


%back substitution  M*x=z
x(n)=z(n)/U(n,n);
for i=n-1:-1:1
   alpha=0;
   for j=i+1:n
      alpha=alpha+x(j)*U(i,j);
   end
   x(i)=(z(i)-alpha)/U(i,i);
end
x=x';
% check
%check=A*x-b
% it should be zero.

% next code is to evaluate inverse of matrix A.
