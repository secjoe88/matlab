% Joey Willhite
% Description:
%     A function to solve systems of first order ODE's using euler's method in
%     two dimensions.
% Inputs:
%     a: Beginning of the interval
%     b: End of interval
%     h: Step size for time interval
%     f1: First ODE in system (as syms function of t, y1, y2) 
%     f2: Second ODE in system (as syms function of t, y1, y2)
%     alpha: Initial condition for first ODE in system
%     beta: Initial condition for second ODE in system
% Outputs:
%     This function has no outputs
%     *Note*
%         Function prints out select values of w_ii's

function euler2d(a, b, h, f1, f2, alpha, beta)
    %set initial values and initial time
    omega=[alpha; beta];
    t=a;
    
    %begin calculating approximations
    for i=1:(b/h)
        %calculate the next approximations in the system
        w_1next=omega(1, i)+h*f1(t, omega(1, i), omega(2, i));
        w_2next=omega(2, i)+h*f2(t, omega(1, i), omega(2, i));
         
        %add new values to omega, selectively print
        omega=double([omega(1,:) w_1next; omega(2, :) w_2next]);
        if mod(t,1)==0
            disp(['i=',int2str(i), '@time=', num2str(t), ': y(t)=', num2str(omega(1, i)),...
                char(9), char(9), 'y''(t)=', num2str(omega(2, i))]);
        end
        
        %increment time
        t=t+h;
        
    end
    
    plot([0:100], omega(1, 1:(1/h):length(omega(1, :))),...
        [0:100], omega(2,1:(1/h):length(omega(2,:))));
end