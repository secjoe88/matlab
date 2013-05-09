% Joey Willhite
% Description:
%     A function to approximate solutions to ODE's using the fourth order 
%     Adams-Bashford method.
% Inputs:
%     a: Beginning of interval of approximation
%     b: End of interval of approximation
%     n: Number of subintervals
%     iv: Vector of initial values w_0, w_1, w_2, and w_3
%     f: Function to approximate (as syms function of t,y)
%     y: Exact solution to ODE (used for calculating error)
% Outputs:
%     

function omega=ABash4(a, b, n, iv, f, y)
    %set value of delta(t) and include initial values in solution
    %approximation
    h=(b-a)/n;
    omega=iv;
    t=a+(3*h);
    
    %calculate succesive approximations w_i at t_i
    for i=4:n
        %calculate w_i+1
        w_next=omega(i) +(h/24)*(55*f(t, omega(i))-59*f(t-h,omega(i-1))+...
            37*f(t-2*h, omega(i-2))-9*f(t-3*h, omega(i-3)));
        %add to list of approximations
        omega=[double(omega) double(w_next)];
        
        %increment t
        t=t+h;
    end
    
    printSolution(a, h, n, omega, y);
end

function printSolution(a, h, n, omega, y)
    format long;
    for i=0:n
        disp(['w_', int2str(i), ' @ t_', int2str(i), '=', num2str(a+i*h), ':', char(9), ...
            num2str(omega(i+1)), char(9), ';error=', num2str(abs(omega(i+1)-double(y(a+i*h))))]);
    end
end