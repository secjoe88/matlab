% Description:
%     A function to approximate solutions to second order nonlinear ODE 
%     y"(x,t)=f(x,y,y') using Newton's Method of non-linear shooting
% Inputs:
%     f:          cell array of symfun {f1,f2;f3,f4} where [f1,f2] is the
%                 2nd order ODE y"=f(x,y,y') written as a first order system
%                 [y1',y2']=[f1(x,y1,y2),f2(x,y1,y2)], f3=df/dy, and f4=df/dy'
%     xbounds:    array [a,b] of boundaries a<=x<=b
%     bcond:      array [alpha, beta] of boundary conditions y(a)=alpha, y(b)=beta
%     dxparam:    parameter for calculating dx. if dxparam<=1, will be treated as a step size
%                 if dx>1 it will be treated as the number of x subintervals
%     tol:        calculation tolerance
%     max_it:     maximum number of allowed iterations
% Output:
%     w:          approximations w(1,i) to y(x_i) and w(2,i) to y'(x_i) for i=1...N
function w=newtonShooting(f,xbounds,bcond,dxparam,tol,max_it)
    a=xbounds(1); b=xbounds(2);alpha=bcond(1);beta=bcond(2);
    %define step size and number of subintervals
    if dxparam<=1
        n=(b-a)/dxparam; dx=dxparam;
    else
        n=dxparam; dx=(b-a)/n;
    end
    index=1; tk=(beta-alpha)/(b-a);w=zeros(2,n+1); u=w;k=zeros(4,2);kp=k;
    
    
    while index<max_it
        %set starting values from initial values
        %u=zeros(size(w));
        w(1,1)=alpha; w(2,1)=tk;u(1,1)=0; u(2,1)=1;
        %begin modified Runge-Kutta method
        for i=2:n+1
            x=a+(i-2)*dx;
            %calculate values of K1-4 equations for system f at point
            %x_(i-1)
            for j=1:2
                k(1,j)=dx*f{1,j}(x, w(1,i-1), w(2,i-1));
            end
            for j=1:2
                k(2,j)=dx*f{1,j}(x+(dx/2), w(1,i-1)+(1/2)*k(1,1), w(2,i-1)+(1/2)*k(1,2));
            end
            for j=1:2
                k(3,j)=dx*f{1,j}(x+(dx/2), w(1,i-1)+ (1/2)*k(2,1), w(2,i-1)+(1/2)*k(2,2));
            end
            for j=1:2
                k(4,j)=dx*f{1,j}(x+dx, w(1,i-1)+k(3,1), w(2,i-1)+k(3,2));
            end
            %use calculated k equations to determine next iteration for
            %y()
            for j=1:2
                w(j,i)=w(j,i-1)+(k(1,j)+2*k(2,j)+2*k(3,j)+k(4,j))/6;
            end
            
            %calculate values of K1-4 equations for the system [z1',z2']=[z2,f3z1+f4z2]
            %used in the newton method
            coeff=[0,1,1,2];
            kp(1,1)=dx*u(2,i-1); kp(1,2)=dx*(f{2,1}(x,w(1,i-1),w(2,i-1))*u(1,i-1)+f{2,2}(x,w(1,i-1),w(2,i-1))*u(2,i-1));
            for l=2:4
                kp(l,1)=dx*(u(2,i-1)+(1/2)*coeff(l)*kp(l-1,2)); 
                kp(l,2)=dx*(f{2,1}(x+(1/2)*dx*coeff(l),w(1,i-1),w(2,i-1))*(u(1,i-1)+(1/2)*coeff(l)*kp(l-1,1))...
                    +f{2,2}(x+(1/2)*dx*coeff(l),w(1,i-1),w(2,i-1))*(u(2,i-1)+(1/2)*coeff(l)*kp(l-1,2)));
            end
            %use calculated kp equations to determine next iteration for
            %z
            u(1,i)=u(1,i-1)+(kp(1,1)+2*kp(2,1)+2*kp(3,1)+kp(4,1))/6;
            u(2,i)=u(2,i-1)+(kp(1,2)+2*kp(2,2)+2*kp(3,2)+kp(4,2))/6;
        end
        %plot(linspace(1,3,33),w(1,:));hold on;plot(linspace(1,3,33),u(1,:));
        %if the approximation is within our tolerance, break the loop
        if abs(w(1,end)-beta)<=tol
            break;
        %otherwise, calculate the next tk using newton's method, and repeat
        else
            tk=tk-(w(1,end)-beta)/u(1,end)
            index=index+1;
        end
    end
    if index>=max_it
        disp('Max iterations exceeded');
    end
    
end