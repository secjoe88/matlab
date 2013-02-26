% Joey Willhite
% Description:
%     simple function to determine the nth order approximation of a derivative 
%     function using Richardson's Extrapolation
% Inputs:
%     x_0: the point being evaluated
%     f: the function being evaluated
%     h: the stepsize
%     iterations: n
% Outputs:
%     result: the nth approximation evaluated at x_0

function result = richExtrap(x_0, f, h, iterations)
    if(iterations ==1),
        result=(f(x_0+h)-f(x_0))/h;
    else
        result = richExtrap(x_0, f, h/2, iterations-1)...
            +((richExtrap(x_0,f, h/2,iterations-1)...
            -richExtrap(x_0, f, h, iterations-1))/((4^(iterations-1))-1));
    end
end
    