% Joey Willhite
% Description:
%     A function to approximate solutions to ODE's using the midpoint method
% Inputs:
%     a: Beginning of interval
%     b: End of interval
%     n: number of subintervals
%     alpha: initial condition
%     f: Function to approximate (as syms function of t, y)
%     y: Exact solution (for error calculation)
% Outputs:
%     omega: Vector of approximations at t_i
%     *Note*
%         Function prints out values of w_i at time t_i
function omega=odemidpoint(a, b, n, alpha, f, y)
    %set step size, fill initial approximation value, set initial time
    h=(b-a)/n;
    omega=alpha;
    t=a;
    
    %begin approximating
    for i=1:n
        %calculate next w value
        w_next=omega(i)+h*f(t+(h/2), omega(i)+(h/2)*f(t, omega(i)));
        %store in list of w values
        omega=[double(omega) double(w_next)];
        %increment t
        t=t+h;
    end
    
    printSolutions(a, h, n, omega, y);
end

function printSolutions(a, h, n, omega, y)
    for i=1:n
        disp(['w_', int2str(i-1), ' @ t_', int2str(i-1), '=', ...
            num2str(a+((i-1)*h)), ':', char(9), num2str(omega(i)), char(9),...
            ';exact=', num2str(double(y(a+((i-1)*h))))]);
    end
end