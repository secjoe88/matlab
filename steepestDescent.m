% Description:
%     A function to approximate solutions to nonlinear systems using the steepest descent method
% Inputs:
%     functions:      A (column) cell array of functions as symfun
%     x0:             A (row) array corresponding to the initial guess
%     n:              Number of equations to solve for
%     terminator:     Terminating case. If 0<terminator<1 it will be treated as a tolerance. If
%                     terminator>1 it will be treated as a max iteration count
% Output:
%     solutions:      Approximation to the solution of the system
function solution=steepestDescent(functions, x0, n, terminator)
    %first, we create our equation g(x1,...,xn) from our inputed functions
    gx=symfun(0,symvar([functions{:}]));
    for i=1:n
        gx=gx+functions{i}^2;
    end
    
    %uncomment this line if you need only to minimize the input
    %gx=functions{1};
    %%begin algorithm
    
    %step3 
    g1=subs(gx,symvar([functions{:}]), x0);
    z=subs(gradient(gx, symvar(gx)), symvar(gx), x0);
    z0=norm(z);
    
    %step4 
    if z0==0
        solution=[x0,g1];
        return
    end
    
    %step 5
    z=z/z0; 
    alpha1=0; 
    alpha3=1; 
    g3=subs(gx,symvar(gx),x0-double(alpha3*z));
    
    %step6
    while (g3>=g1)
        %step 7
        alpha3=alpha3/2;
        g3=subs(gx,symvar(gx),x0-double(alpha3*z));
        %step8 (serves also as terminating condition for tolerance terminator)
        if terminator<1
            if alpha3<terminator/2
                solutions=[x0 g1];
            end
        end
    end
    %step9
    alpha2=alpha3/2; g2=subs(gx,symvar(gx),x0-double(alpha2*z));
    %step10
    h1=(g2-g1)/alpha2; h2=(g3-g2)/(alpha3-alpha2); h3=(h2-h1)/alpha3;
    %step11
    alpha0=0.5*(alpha2-h1/h3); g0=subs(gx, symvar(gx),x0-double(alpha0*z));
    %step12/13
    if g0<g3
        x0=x0-double(alpha0*z)
    else
        x0=x0-double(alpha3*z)
    end
    %step14
    g=min([double(g0), double(g3)]);
    if terminator<1 
        if norm(g-g1)<terminator
            solution=[x0, g]
            return
        else
            %if our next guess is not within tolerance, recurse
            steepestDescent(functions, x0, n, terminator)
        end
    else
       %%terminating condition for iterative terminator
        terminator=terminator-1;
        if terminator==0
            solutions=[x0, g]
        else
            steepestDescent(functions, x0, n, terminator)
        end
    end
    
        
    