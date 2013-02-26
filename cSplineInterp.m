% Joey Willhite
% Inputs:
%     x=array of x values [x_0,...,x_n]
%     a=array of a=f(x_n) values [a_0,..,a_n]
% This script has no outputs
% *Note*
%     Script prints rows a_j, b_j, c_j, d_j

function cSplineInterp(x, a)
    n=size(x);
    n=n(2);
    
    %create variables (in order of appearance)
    h=zeros(1, n);
    alpha=zeros(1, n-1);
    l=zeros(1, n);
    mu=zeros(1, n-1);
    z=zeros(1,n);
    c=zeros(1,n);
    b=zeros(1,n-1);
    d=zeros(1,n-1);
    
    %set h and alpha values (Step 1/2)
    h(1)=x(2)-x(1);
    for i=2:n-1,
        h(i)=x(i+1)-x(i);
        alpha(i)=(3/h(i))*(a(i+1)-a(i))-(3/h(i-1))*(a(i)-a(i-1));
    end
    
    %Step 3 set initial l, mu, z values
    %(initial mu and z values already set)
    l(1)=1;
    
    %Step 4
    for i=2:n-1,
        l(i)=2*(x(i+1)-x(i-1))-h(i-1)*mu(i-1);
        mu(i)=h(i)/l(i);
        z(i)=(alpha(i)-h(i-1)*z(i-1))/l(i);
    end
    
    %Step 5 set l,z,c values at n
    l(n)=1;
    z(n)=0;
    c(n)=0;
    
    %Step 6 calulate c_j, b_j, d_j
    for i=n-1:-1:1,
        c(i)=z(i)-mu(i)*c(i+1);
        b(i)=(a(i+1)-a(i))/h(i)-h(i)*(c(i+1)+2*c(i))/3;
        d(i)=(c(i+1)-c(i))/(3*h(i));
    end
    
    disp(a);
    disp(b);
    disp(c);
    disp(d);
    
    
    
end