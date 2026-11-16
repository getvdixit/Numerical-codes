% Romberg Integration
%define a function f. ep is the tolerance, m is the number of iterations, b and a are 
% the upper and lower limits respectively
function a =Romberg ( ep ,m ,b , a)
f=inline('x.^4');
%f=inline('(1/(exp(1)-1))*exp(x)');
r=zeros(m);
k=1;h=b-a;		%initialize the variables
r(1,1)=h*(f(a)+f(b))/2;	
while(k<m)
    k=k+1;           h=h/2;	
    c=sum(f(a+h:2*h:b-h));
    
    r(k,1)=0.5*(r(k-1,1)+2*h*c)			%use the trapezoidal rule
                              			%hint: you can generate the vector a:h:b and use sum
    for j=2:k
    r(k,j)=(4^(j-1)*r(k,j-1)-r(k-1,j-1))/(4^(j-1)-1)			%use extrapolated estimates 
    end
    e=r(k,k)-r(k-1,k-1);			% estimate truncation error
    if( abs(e) < ep)
        a=r(k,k);
        return			%break out of the while loop and return to command window 
    end
end
a=r(k,k);				%use format command for displaying more digits after decimal

%Assignment:	You can test the above program for the function f(x)= c exp(x)  
%for c=1/(e-1)  , a=0 and b=1 the answer should converge to 1.
