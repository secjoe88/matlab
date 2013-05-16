% Joey Willhite
% Description:
%     A function that uses the linear shooting method in combination with the 
%     Runge Kutta fourth order method for approximating solutions to second order
%     boundary value problems.
% Inputs:
%     a: Beginning of interval
%     b: End of interval
%     alpha: Beginning boundary conditions
%     beta: End boundary condition
%     n: number of subintervals
%     f: Function list as 2x2 cell array. 
%     
%         Function list is formatted with f_1i being the first system used 
%         by the linear shooting method derived by rewriting the original second order DE 
%         (including the r(t) term) as a system of 2 first order DE's, 
%         and f_2i being the second system used by the linear shooting method 
%         derived by removing the r(t) term and rewriting the second order
%         DE as a system of 2 first order DE's.
%         
% Outputs:
%     omegas: Approximations omega(1, i) to y(t) and omega(2, i) to y'(t) for
%             t in [a, b].
%     *Note* 
%     Function also prints out select values of Runge Kutta evaluations, select
%     values of calculated approximation to second order DE.

function linearShooting(a, b, alpha, beta, n, f, exactsol)
    %set initial value for h
    h=(b-a)/n;
    
    %fill u, v, with Runge Kutta approximations 
    disp(['RK4 Approximation for y"=p(t)y''+q(t)y+r(t)']);
    u=RK4systems(a, b, 2, h, f(1,:),[alpha 0]);
    disp('RK4 Approximation for y"=p(t)y''+q(t)y');
    v=RK4systems(a, b, 2, h, f(2,:),  [0 1]);
    
    %use RungeKutta approximations to calculate linear shooting
    %approximnations
    omega=zeros(2, n+1);
    omega(1,1)=alpha; 
    omega(2,1)=(beta-u(1,n+1))/v(1,n+1);
    
    disp(['Combined RK4 Approximations with linear shooting', char(13), ...
        'i=1@time t0=', num2str(a) ,': y(t)=', num2str(omega(1,1)), char(9), 'y''(t)=', ...
        num2str(omega(2,1))]);
    for i=2:n+1
        omega(1,i)=u(1,i)+omega(2,1)*v(1, i);
        omega(2,i)=u(2,i)+omega(2,1)*v(2,i);
        
        if mod(i-1,2)==0
            resultmessage=['i=', num2str(i), '@time t',num2str(i-1), '=', ...
                num2str((i-1)*h) ,': ~y(t)=', num2str(omega(1,i)), char(9), ...
                '~y''(t)=', num2str(omega(2,i))];
            if nargin==7
                resultmessage=[resultmessage, char(13), char(9),'y(t)=', ...
                    num2str(double(exactsol((i-1)*h))),';', char(9), ' E(y(t))='...
                    , num2str(double(abs(omega(1,i)-exactsol((i-1)*h))))];
            end
            disp(resultmessage);
        end
        
    end
    plot([a:h:b], omega(1,:), [a:h:b], omega(2,:));
end