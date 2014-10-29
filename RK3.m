% Joey Willhite
% Description:
%     A function to approximate solutions to ODE's using the Runge-Kutta
%     method of order 3
% Inputs:
%     a: Beginning of interval
%     b: End of interval
%     n: Number of sub-intervals
%     alpha: Initial condition
%     f: Function to approximate
% Outputs:
%     omega: List of approximations for t_i
%     *Note*
%         Function prints out approximation and exact value at t_i
function omega=RK3(a, b, n, alpha, f, y)
    %set step size, fill in initial value, set t_0
    h=(b-a)/n;
    omega=alpha;
    t=a;
    
    %begin calculating approximations
    for i=1:n
        %calculate next approximation using RK3
        w_next=omega(i)+(h/4)*(f(t, omega(i))+3*f(t+(2*h)/3, omega(i)+...
            ((2*h)/3)*f(t+(h/3), omega(i)+(h/3)*f(t, omega(i)))));
        %add next approximation to list of approximations
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