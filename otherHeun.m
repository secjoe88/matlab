% Joey Willhite
% Description:
%     A function to solve first order linear ODEs using the Taylor method of
%     order 
% Inputs:
%     diffeq:     symfun f(t,y) of the ODE y'=f(t,y)
%     domain:     an array [a,b] of the domain of approximation
%     ya:         initial value y(a)=ya
%     h:          step size
%                 *NOTE*  if h<=1, it will be treated as a step size, 
%                 otherwise it will be treated as a mesh point count
% Outputs:
%     solutions:  approximation to solution
function solution=otherHeun(diffeq, domain, ya, h)
    %code to make function robust to step sizes and mesh points
    if h>1
        %treat h like a mesh point count
        index=h;
        h=(domain(2)-domain(1))/h;
    else
        %treat h like a step size
        index=(domain(2)-domain(1))/h;
    end
    %setup variables for Heun method
    t=domain(1);
    w=ya;
    solution=[t,w];
    
    %%commence Heun method
    for i=1:index
        q1=subs(diffeq,symvar(diffeq), solution(i,1:length(symvar(diffeq))));
        q2=diffeq(t+2*h/3, w+2*h/3*q1);
        
        w=double(w+(h/4)*(q1+3*q2));
        t=domain(1)+i*h;
        solution=[solution; t,w];
    end
end