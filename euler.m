% Joey Willhite
% Description:   
%     Function to approximate solutions to an ordinary differential equation
%     on a finite interval using Euler method. Function also
%     prints exact solution based on inputed particular solution
% Inputs:
%     a: beginning of interval
%     b: end of interval
%     N: number of sub-intervals
%     alpha: initial condition y(a)=alpha
%     f: ODE to approximate (as syms)
%     y: particular solution to ODE
% Outputs:
%     Function has no outputs
%     *Note* Function prints approximation and exact solution at subsequent
%     values of t


function euler(a,b,N, alpha, f, y)
    h=(b-a)/N;
    t=a;
    w=alpha;
    disp(['Approximation w_0=', num2str(alpha)]);
    
    for i=1:N
        w=w+(h*eval(f(t,w)));
        t=a+(i*h);
        disp(['Approximation w_', int2str(i), ' at t_', int2str(i), '=',...
            num2str(t),': ', num2str(w), ',', 9,9, 'exact y(t_', int2str(i), ...
            ')= ', num2str(eval(y(t)))]);
    end
end