% Joey Willhite
% Description:
%     A function to approximate solutions to ODE's using the second order 
%     Adams-Bashford method.
% Inputs:
%     a: Beginning of interval of approximation
%     b: End of interval of approximation
%     n: Number of subintervals
%     iv: Vector of initial values w_0, and w_1,
%     f: Function to approximate (as syms function of t,y)
%     y: Exact solution to ODE (used for calculating error)
% Outputs:
%     

function ws=ABash2(a, b, n, iv, f, y)
    %set value of delta(t) and include initial values in solution
    %approximation
    h=(b-a)/n;
    ws=iv;
    t=a+(3*h);
    
    %calculate succesive approximations w_i at t_i
    for i=2:n
        %calculate w_i+1
        w_next=ws(i) +(h/2)*(3*f(t, ws(i))-f(t-h,ws(i-1)));
        %add to list of approximations
        ws=[double(ws), double(w_next)];
        
        %increment t
        t=t+h;
    end
    
    printSolution(a, h, n, ws, y);
end

function printSolution(a, h, n, omega, y)
    format long;
    for i=0:n
        disp(['w_', int2str(i), ' @ t_', int2str(i), '=', num2str(a+i*h), ':', char(9), ...
            num2str(omega(i+1)), char(9), ';error=', num2str(abs(omega(i+1)-double(y(a+i*h))))]);
    end
end