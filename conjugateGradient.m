% Joey Willhite
% Description:
%     A method to compute the solution to a linear system using the preconditioned
%     conjugate gradient method to within a specified tolerance or specified
%     number of iterations.
% Inputs:
%     n: Number of unknowns/equation in system
%     A: Coefficient matrix
%     b: Vector of column n+1 or augmented coefficient matrix. 
%         *Note* must be entered as a column vector
%     C: Preconditioning matrix
%     x: Initial vector approximation.
%         *Note* must be entered as a column vector
%     max_iterations: Desired maximum iterations, if 0 continue until
%       tolerance is reached
%     tol: Desired approximation tolerance
% Outputs:
%     This function has no outputs.
%     *Note* Function prints out resulting approximation along with resulting
%         residual vector.
%	also, git test
function x= conjugateGradient(n, A, b, C, x, max_iterations, tol)
    %Step 1: Establish initial residual; residual with preconditioning;
    %instantiate necessary variables
    r=b-(A*x); 
    w=C*r;
    v=transpose(C)*w;
    alpha=transpose(w)*w;
    u=zeros(n,1);
    beta=0;
    s=0;
    
    %Step 2-3: Begin iteration
    k=1;
    while(k<=max_iterations || ~max_iterations)
        %Step 4: Check if initial value is within tolerance, if so return
        if(inf_norm(v, n)<=tol)
            disp(['Solution vector:[' num2str(transpose(x)) ']']);
            disp(['Residual vector:[' num2str(transpose(r)) ']']);
            disp(['Iterations: ' num2str(k)]);
            return;
        end
        
        %Step 5: Compute next iteration
        u=A*v;
        t=alpha/(transpose(u)*v);
        x=x+(t*v);
        r=r-(t*u);
        w=C*r;
        beta=transpose(w)*w;
        
        %Step 6: Check if next approximation is within tolerance, if so
        %return, otherwise resume with calculations
        if(abs(beta)<=tol)
            if(inf_norm(r, n)<=tol)
            disp(['Solution vector:[' num2str(transpose(x)) ']']);
            disp(['Residual vector:[' num2str(transpose(r)) ']']);
            disp(['Iterations: ' num2str(k)]);
            return;
            end
        end
        
        %Step 7: Prepare for next iteration
        s=beta/alpha;
        v=(transpose(C)*w)+(s*v);
        alpha=beta;
        k=k+1;
    end
    disp('Maximum iterations exceeded');
end

function norm=inf_norm(v, n)
    norm=0;
    for i=1:n
        if(abs(v(i))>=norm)
            norm=abs(v(i));
        end
    end
end