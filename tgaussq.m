%ASS==Please compute all the non-zero mass moment of inertia of 
%(a) a right angled tetrahedron with sides a=1(say,x),b=1(say,y),c=2,d=3(height i.e. z).
%(b) a hemisphere of radious 2 by gauss quadrature with 2,3,4 and 5 gauss points in each direction.
% compare answers with standard analytical expression.
%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function integrates a 3-d integral for n-point gauss-quadrature. With  %%
%% limits x=(a , b) ; y=(g1(x) , g2(x)) ; z=(h1(x,y) , h2(x,y))  for the funct.%%
%% f(x,y,z)dxdydz. f=function not taking transformations into account.         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tgaussq(n)
x0=1; y0=1; z0=3;
a=-1;  b=1;
g1=inline('0','x'); g2=inline('1-x','x');              %limits g2(x)=y0*(x0-x)/x0
h1=inline('0','x','y'); h2=inline('3*(x+y)','x','y');      %limits h2(x,y)=z0*(x0*y+y0*x)/(x0*y0)

f=inline('exp(x)*sin(x)/sqrt(1-x^2)','x','y','z'); % function to be integrated.
switch n      %the array of evaluation points (Use a standard table)
case 1
   w=[0];
case 2
   w=[-0.5773502691896257  0.5773502691896257];
case 3
   w=[-0.77459666924148  0  0.77459666924148];
case 4
   w=[-0.8611363115940   -0.339981043584   0.339981043584   0.8611363115940];
case 5
   w=[-0.9061798  -0.538469310   0     0.538469310    0.9061798];
end

switch n               %the array of coefficients
case 1
   c=[2];
case 2
   c=[1 1];
case 3
   c=[5/9 8/9 5/9];
case 4
   c=[0.34785 0.65214  0.65214 0.34785];
case 5
   c=[0.2369  0.4786  0.5688  0.47862   0.2369];
end

for itr = 1 : n  %scaling the function
   x(itr) = ( ( 1.0 - w(itr) ) * a + ( 1.0 + w(itr) ) * b ) / 2.0;
end

 %fvec=f(x);
%a=(fvec*c'/2)*(b-a)		%the approximate integral


sum=0;                  % INTEGRAL FOR MOMENT OF INERTIA.
for i=1:n
   for itr = 1 : n        % scaling the points.
      y(itr) = ( ( 1.0 - w(itr) ) * g1(x(i)) + ( 1.0 + w(itr) ) * g2(x(i)) ) / 2.0;	 
   end
   
   for j=1:n
      for itr = 1 : n                %scaling the points
         z(itr) = ( ( 1.0 - w(itr) ) * h1(x(i),y(j)) + ( 1.0 + w(itr) ) * h2(x(i),y(j)) ) / 2.0;	 
      end
        
      for k=1:n
         mult=(b-a)*(g2(x(i))-g1(x(i)))*(h2(x(i),y(j))-h1(x(i),y(j)))/8;
         sum=sum+c(i)*c(j)*c(k)*f(x(i),y(j),z(k))*mult;
      end
   end
end

disp('The moment of inertia is: ')
disp(sum);


%%%%%%%%%%%%%%%NOTE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For integral less than three degree use other limits like g,h   %%
%% from 0 to 1. keep function as it is you will get correct answer.%%
%% In case of constant limits keep function constant.              %%
%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         