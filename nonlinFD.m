% Description:
%     A function to approximate solutions to the nonlinear boundary value ODE 
%     u"=f(x,u,u') using the Newton's iterative finite difference method
% Inputs:
%     f:          cell array of symfun {f1,f2,f3} where for the BVP u"=f(x,u,u'),
%                 f1=f, f2=df/dy, and f3=df/dy'
%     xbounds:    array [a,b] of boundaries a<=x<=b
%     bconds:     array [alpha,beta] of boundary conditions u(a)=alpha, u(b)=beta
%     N:          number of subintervals
%     tol:        tolerance
%     max_it:     maximum iterations
% Output:
%     w:          array of approximations w(i) to solutions y(x_i) for i=1...N
function w=nonlinFD(f,xbounds,bcond,N,tol,max_it)
    a=xbounds(1); b=xbounds(2);alpha=bcond(1);beta=bcond(2);w=zeros(1,N+2);
    dx=(b-a)/(N+1);
    w(1)=alpha;w(N+2)=beta;
    for i=2:(N+1)
        w(i)=alpha+(i-1)*((beta-alpha)/(b-a))*dx;
    end
    k=1;
    ai=zeros(1,N);bi=ai;ci=ai;di=ai;
    while k<=max_it
        %step 5
        x=a+dx;
        t=(w(3)-alpha)/(2*dx);
        ai(1)=2+(dx^2)*f{2}(x,w(2),t);
        bi(1)=-1+(dx/2)*f{3}(x,w(2),t);
        di(1)=-(2*w(2)-w(3)-alpha+(dx^2)*f{1}(x,w(2),t));
        %step6
        for i=2:(N-1)
            x=a+i*dx;
            t=(w(i+2)-w(i))/(2*dx);
            ai(i)=2+(dx^2)*f{2}(x,w(i+1),t);
            bi(i)=-1+(dx/2)*f{3}(x,w(i+1),t);
            ci(i)=-1-(dx/2)*f{3}(x,w(i+1),t);
            di(i)=-(2*w(i+1)-w(i+2)-w(i)+(dx^2)*f{1}(x,w(i+1),t));
        end
        %step 7
        x=b-dx;
        t=(beta-w(N))/(2*dx);
        ai(N)=2+(dx^2)*f{2}(x,w(N+1),t);
        ci(N)=-1-(dx/2)*f{3}(x,w(N+1),t);
        di(N)=-(2*w(N+1)-w(N)-beta+(dx^2)*f{1}(x,w(N+1),t));
        
        %step 8 (solve tridiagonal system)
        l=ai(1);
        u(1)=bi(1)/ai(1);
        z(1)=di(1)/l(1);
        %step 9
        for i=2:(N-1)
            l(i)=ai(i)-ci(i)*u(i-1);
            u(i)=bi(i)/l(i);
            z(i)=(di(i)-ci(i)*z(i-1))/l(i);
        end
        %step 10
        l(N)=ai(N)-ci(N)*u(N-1);
        z(N)=(di(N)-ci(N)*z(N-1))/l(N);
        
        %step 11
        v=zeros(1,N); v(N)=z(N);
        w(N+1)=w(N+1)+v(N);
        %step 12
        for i=(N-1):-1:1
            v(i)=z(i)-u(i)*v(i+1);
            w(i+1)=w(i+1)+v(i);
        end
        %step 13
        if (norm(v,2)<=tol)
            for i=1:length(w(1,:))
                w(2,i)=a+(i-1)*dx;
            end
            break;
        end
        k=k+1;
        
    end
end

    
    
    
    
    
    
    
    
    
    