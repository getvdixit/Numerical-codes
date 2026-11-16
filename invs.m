% This function finds the inverse of A using Gauss method (gauss.m).

function B=invs(A)
[m,n]=size(A);
b=eye(n);
for i=1:n
   b(:,i)=gauss(A,b(:,i));
end
B=b;
    
   