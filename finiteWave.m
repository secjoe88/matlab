% Description:
%     A function to solve wave equations using the wave equation finite difference
%     algorithm
% Inputs:
%     l:      endpoint 0<x<l
%     T:      maximum time 0<t<T
%     alpha:  constant in the equation u_tt-(alpha^2)u_xx=0
%     h:      x step size
%     k:      t step size
%     f:      initial condition u(x,0)=f(x) as symfun
%     g:      initial condition u_t(x,0)=g(x) as symfun
% Outputs:
%     xs:     array of x values x0+ih for i=0,...,m
%     ys:     array of t values t0+jk for j=0,...,n
%     w:      mxn matrix of values x_ij approximations to solution of pde
function [xs,ts,w]=finiteWave(l,T,alpha,h,k,f,g)
    %Step 1
    m=floor(l/h);
    n=floor(T/k);
    lambda=(k*alpha)/h;
    w=zeros(m+1,n+1);
    %Step2
    for j=1:n
        w(1,j+1)=0;
        w(m+1,j+1)=0;
    end
    %Step 3
    w(1,1)=f(0);
    w(m+1,1)=f(l);
    %Step 4
    for i=1:m-1
        w(i+1,1)=f(i*h);
        w(i+1,2)=(1-lambda^2)*f(i*h)+((lambda^2)/2)*(f((i+1)*h)+f((i-1)*h))+k*g(i*h);
    end
    %Step 5
    for j=1:(n-1)
        for i=1:(m-1)
            w(i+1,j+2)=2*(1-lambda^2)*w(i+1,j+1)+(lambda^2)*(w(i+2,j+1)+w(i, j+1))-w(i+1,j);
        end
    end
    %Step 6
    xs=0;
    ts=0;
    for i=1:m
        xs=[xs, xs(i)+h];
    end
    for i=1:n
        ts=[ts; ts(i)+k];
    end
        
end