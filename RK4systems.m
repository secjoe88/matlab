% Description:
%     A function to approximate solutions to 2nd order systems of ODE's using
%     Runge Kutta Method in 2-dimensions
% Inputs:
%     a: Beginning of interval
%     b: End of interval
%     m: Number of equations
%     h: Step size
%     f: Cell array of equations of system
%     alphas: Vectors of initial values for equations
% Outputs:
%     This function has no outputs
%     *Note*
%         Function prints out values of eq1 and eq2 at various values of t and
%         plots both equations as functions of t

function omegas=RK4systems(a, b, m, h, f, alphas)
    %set initial time, create k's, create w's
    t=a;
    k=zeros(4, 2);
    w=zeros(m,1);
    omegas=zeros(m,0);
    %fill initial values 
    for j=1:m
        w(j)=alphas(j);
        omegas(j,1)=alphas(j);
    end
    
    %begin calculating approximations
    for i=1:(b/h)
        %calculate values of K1-4 equations
        for j=1:m
            k(1,j)=h*f{j}(t, w(1), w(2));
        end
        for j=1:m
            k(2,j)=h*f{j}(t+(h/2), w(1)+(1/2)*k(1,1), w(2)+(1/2)*k(1,2));
        end
        for j=1:m
            k(3,j)=h*f{j}(t+(h/2), w(1)+ (1/2)*k(2,1), w(2)+(1/2)*k(2,2));
        end
        for j=1:m
            k(4,j)=h*f{j}(t+h, w(1)+k(3,1), w(2)+k(3,2));
        end
        
        %use K equations to calculate next approximation to the equations
        %in the systems
        for j=1:m
            w(j)=w(j)+(k(1,j)+2*k(2,j)+2*k(3,j)+k(4,j))/6;
        end
        %store next calculated approximation
        omegas=double([omegas w]);
        %print selectively
        if mod(t,1)==0
            disp(['i=',int2str(i), '@time=', num2str(t), ': y(t)=', num2str(omegas(1, i)),...
                char(9), char(9), 'y''(t)=', num2str(omegas(2, i))]);
        end
        t=t+h;
    end
    
    plot([0:100], omegas(1, 1:(1/h):length(omegas(1, :))),...
        [0:100], omegas(2,1:(1/h):length(omegas(2,:))));
end
