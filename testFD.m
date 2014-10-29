% Summary:
%     A function to approximate solutions to the PDE u_t=-uu_x+eu_xx+f(x,t) for 
%     e=1 and f(x,t)=sin(pi*x)*(pi*cos(t)^2*cos(pi*x)+e*pi^2*cos(t)-sin(t))
%     and initial condition u(x,t)=cos(t)*sin(pi*x) on the bounds
%     -1<x<1,0<t<2. It employs a 3-tier semi-implicit finite difference
%     scheme, and uses gauss-seidel to solve the resulting implicit system.
%     Used for problem #3 of the final to test convergence of the scheme.
% Inputs: 
%     subs:   the temporal/spatial step size (for testing purposes, dx=dt)
%     tol:    approximation tolerance (stopping condition)
%     it:     maximum number of GS iterations (if 0, will be treated as infinite)
% Outputs:
%     meshx:      meshgrid of x values used for plotting. (created by the meshgrid function)
%     mesht:      meshgrid of t values used for plotting. (created by the meshgrid function)
%     w:          (p+1)x(m+1) matrix of approximations w(i,j) to the solutions at points x_ij
%                 for i=1,..,(p+1), j=1,...,(m+1)
%     frames:     array of frames of the plot of the error at each iteration
%                 of GS, (for error analysis and debugging purposes)
function [meshx,mesht,w,frames]=testFD(subs,tol,it)
    %variables (bounds,subintervals,etc.)
    a=-1;b=1;c=0;d=2;dt=subs;dx=dt/4;e=1;p=floor((b-a)/dx);m=floor((d-c)/dt);
    e=1;L=(2*dt)/(dx^2);G=dt/dx;beta=e*L;mu=3+(2*beta);
    x=transpose([a:dx:b]);t=transpose([c:dt:d]);
    
    %for testing purposes we set up a grid matrix with exact values along
    %the borders and 0 at all the interior points. we also calculate the
    %exact solution in entirety for error analysis purposes
    [meshx,mesht]=meshgrid(x,t);frames=[];
    w=zeros(p+1,m+1);exact=cos(mesht).*sin(pi*meshx);
    w(:,1)=cos(0).*sin(pi*x);w(:,m+1)=cos(2)*sin(pi*x);
        
    %uncomment to use hardware acceleration when rendering the plot of the
    %error
    %set(gcf,'Renderer', 'opengl');
        
    %in the interest of simply testing the method, we make our initial condition
    %(w:,2) values from the exact solution
    w(:,2)=cos(t(2)).*sin(pi*x);
    
    %now we execute our finite difference scheme on the environment we set
    %up, using gauss seidel to solve the system, and test for convergence
    l=1;
    while true
        for j=2:m-1
            for i=2:p
                %calculate values from boundary conditions and previous time steps
                Cx=-w(i,j-1)+G*w(i,j-1)*(w(i+1,j-1)-w(i-1,j-1))-2*dt*f(x(i),t(j-1));
                Bx=4*w(i,j)-2*G*w(i,j)*(w(i+1,j)-w(i-1,j))+4*dt*f(x(i),t(j));
                ax=beta*(w(i+1,j+1)+w(i-1,j+1));
                %calculate next approximation using gauss-seidel
                z=(1/mu)*(ax+Bx+Cx);
                if i==2 && j==2
                    norm=abs(z-w(i,j+1));
                end
                %residual vector analysis
                if norm<abs(z-w(i,j+1))
                    norm=abs(z-w(i,j+1));
                end;
                %store next approximation
                w(i,j+1)=z;
                %surf(meshx,mesht,w');
            end
        end
        %plots the accuracy of the current iteration as a 3-d surface , and 
        %saves a frame of it to be viewed in a time lapse movie of the error
        %evolution per iteration (for error analysis and debugging purposes)
%         surf(meshx,mesht,w');
%         axis([-1,1,0,2,-1,1]);
%         frames=[frames;getframe()];
        
        %print current norm and current gauss-seidel iteration count
        if l>1
            fprintf(repmat(sprintf('\b'),1,mlen));
        end
        msg=sprintf('l=%d\tnorm=%d',l,norm);mlen=length(msg+1);
        fprintf(msg);
        %check norm against tolerance and iteration count (terminating
        %conditions)
        if norm<=tol || l==it
            disp(' ');%fprintf(repmat(sprintf('\b'),1,mlen));
            break;
        end
        %increment iteration count
        l=l+1;
    end
        
    %helper function outlining the forcing function f(x,t) of the PDE F(u)=f(x,t). 
    %Must be coded explicitly for the purposes of minimizing computational
    %overhead
    function sol=f(x,t)
        sol=sin(pi*x)*(pi*cos(t)^2*cos(pi*x)+e*pi^2*cos(t)-sin(t));
    end
    
    
end