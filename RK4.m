% Joey Willhite
% Description:
%     Function to approximate solutions to an ordinary differential equation
%      on a finite interval using Rung-Kutta method of order 4. Function also
%      prints exact solution based on inputed particular solution
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

function omega=RK4(a, b, N, alpha, f, y)
    h=(b-a)/N;
    t=a;
    w=alpha;
    omega=alpha;
    disp(['Approximation w_0=', int2str(alpha)]);
    k1=0; k2=0; k3=0; k4=0; 
    
    for i=1:N
        k1=h*f(t, w);
        k2=h*f(t+(h/2), w+(k1/2));
        k3=h*f(t+(h/2), w+(k2/2));
        k4=h*f(t+h, w+k3);
        
        w=eval(w+(k1+2*k2+2*k3+k4)/6);
        omega=[double(omega) double(w)];
        t=a+i*h;
        disp(['Approximation w_', int2str(i), ' at t_', int2str(i), '=',...
            num2str(t),': ', num2str(w), ',', 9, 'exact y(t_', int2str(i), ...
            ')= ', num2str(eval(y(t)))]);
    end
end