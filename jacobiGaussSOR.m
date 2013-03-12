% Joey Willhite
% Description:
%     A function to aproximate solutions to a linear system of equations within
%     a specified tolerance or given a maximum number of iterations using 
%     Jacobi's, Gauss-Seidel, or SOR iterative methods.
% Inputs:
%     A: the nxn coefficient matrix of the unknonws of the linear system to be 
%         approximated
%     b: the (n+1) column of the system to be approximated
%     x_0: an initial approximation to the unknowns of the linear system
%     n: the number of unknowns in the system to be approximated
%     tol: the tolerated error of the approximation
%     max_iterations: the maximum alloted iterations for the algorithm to run
%         before returning an approximation. If '0', algorithm will run
%         until approximation matches tolerance
%     option: an integer parameter. if 0, it will calculate using the Jacobi
%         method. if 1 it will calculate using the Gauss-Seidel method. if
%         2 it will calculate using the SOR method.
% Outputs:
%     x: An approximation to the unknowns of the linear system
function x=jacobiGaussSOR(A, b, x_0, n, tol, max_iterations, option)
    k=1;
    x=zeros(1, n);
    if (option==2)
        omega=input('Please input omega value');
    end
    while(k<=max_iterations || ~max_iterations)
        %Step 3, set values of x approximation
        for i=1:n
            switch option
                case 0
                    %Use Jacobi algorithm
                    x(i)=(1/A(i, i))*(-sigma(i, 1, n, A, x_0)+b(i));
                case 1
                    %Use Gauss-Seidel Algorithm
                    x(i)=(1/A(i, i))*(-sigma(i, 1, i-1, A, x)...
                        -sigma(i, i+1, n, A, x_0)+b(i));
                case 2
                    %Use SOR algorithm
                    x(i)=(1-omega)*x_0(i)+(1/A(i, i))*(omega*(-sigma(i, 1, i-1, A, x)...
                        -sigma(i, i+1, n, A, x_0)+b(i)));
                    
                otherwise
                    disp('no valid option');
            end
        end
        %Step 4: computer the error between subsequent iterations. if the 
        %error is within the tolerance, break the operation and return the
        %estimate
        if computeError(x, x_0, n)<=tol
            break;
        end
        %Step 5: increment max iteration counter
        k=k+1;
        %Step 6: set previous estimate as current estimate and repeat
        for i=1:n
            x_0(i)=x(i);
        end
    end
    if(k>=max_iterations && max_iterations)
        disp('Max iterations exceeded');
    end
end

%A function to sum the terms necessary to calculate the sigma in step 3
function result =sigma(row, j, n, A, x)
    result=0;
    for i=j:n
        if i~=row
            result=result+(A(row,i)*x(i));
        end
    end
end
% A function to use the infinity norm to calculate the error between 
% subsequent iterations of the extimated result
function error=computeError(x_1, x_2, n)
    error=0;
    for i=1:n
        if(abs(x_1(i)-x_2(i))>error)
            error=abs(x_1(i)-x_2(i));
        end
    end
end