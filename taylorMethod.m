% Joey Willhite
% Description:
%     A function to solve first order linear ODEs using the Taylor method of
%     order n
% Inputs:
%     diffeq:     symfun f(t,y) of the ODE y'=f(t,y)
%     domain:     an array [a,b] of the domain of approximation
%     ya:         initial value y(a)=ya
%     h:          step size
%                 *NOTE*  if h<=1, it will be treated as a step size, 
%                 otherwise it will be treated as a mesh point count
%     n:          Taylor method order
% Outputs:
%     solutions:  approximation to solution
function apprxs=taylorMethod(diffeq, domain, ya, h, n)
    format long;
    
    %code to make function robust to step sizes and mesh points
    if h>1
        %treat h like a mesh point count
        index=h;
        h=(domain(2)-domain(1))/h;
    else
        %treat h like a step size
        index=(domain(2)-domain(1))/h;
    end
    %prepare variables for approximation
    t=domain(1);
    w=ya;
    apprxs=[t,w];
    
    %generate n order taylor polynomial
    T=taylorPoly(diffeq,n, h);
    
    %create approximations with taylor's method
    for i=1:index
        w=w+double(h*subs(T, symvar(T), apprxs(i,1:length(symvar(T)))));
        t=domain(1)+i*h;
        apprxs=double([[apprxs]; [t,w]]);
    end
end


%helper function for creating taylor polynomials
function result=taylorPoly(diffeq, order, h)
    vars=symvar(diffeq);
    temp=diffeq;
    result=diffeq;
    for i=1:order-1
        %generate ith taylor term
        temp=diff(temp,vars(1))+diffeq*diff(temp,vars(2));
        %add it to result
        result=result+((h^i)/factorial(i+1))*temp;
    end
    
end