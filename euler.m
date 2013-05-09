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
%     solution: Vector containing approximations at t_i's
%     *Note* Function prints approximation and exact solution at subsequent
%     values of t


function solution=euler(a,b,N, alpha, f, y)
    h=(b-a)/N;
    t=a;
    w=alpha;
    disp(['Approximation w_0 at t_0=', num2str(a), char(9), num2str(a) ]);
    solution=alpha;
    
    format long; 
    for i=1:N
        w=w+(h*eval(f(t,w)));
        t=a+(i*h);
        disp(['Approximation w_', int2str(i), ' at t_', int2str(i), '=',...
            num2str(t),': ', char(9), num2str(w), ';', 9,9, 'exact y(t_', int2str(i), ...
            ')= ', num2str(eval(y(t)))]);
        solution=[solution w];
    end
end