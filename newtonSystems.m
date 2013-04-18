% Joey Willhite
% Description:
%     A function to solve nonlinear systems using the Newton method in 2 dimensions
% Inputs:
%     func: System of equations expresed as a vector of symbolic expressions
%     init: Initial approximation as a column vector
%     tol: Tolerance
%     nMax: Max iterations (input -1 to iterate indefinitely)
% Outputs:
%     x: Approximation to system solution
function x=newtonSystems(func, init, tol, nMax)
    syms x1 x2;
    %create inline expressions for F and Jacobian of F
    F=inline(func);
    J=inline(jacobian(func));
    %check if max iterations are exceeded, if so return a message
    %indicating such
    if nMax==0
        disp('Max iterations exceeded');
    end
    
    %if the current approximation is a solution, return it; otherwise
    %calculate the next iteration
    if F(init(1), init(2))==0
        x=init;
        return;
    else
        y=linsolve(J(init(1), init(2)), -F(init(1), init(2)));
        next=init+y;
    end
    
    %if the next iteration is within tolerance, return it. otherwise,
    %recurse
    if(norm((next-init), Inf)<=tol)
        x=next;
        return;
    else
        newtonSystems(func, next, tol, nMax-1)
    end
end