% Summary:
%     A function to approximate solutions to the PDE u_t=-uu_x+eu_xx for 
%     e=.1 and and initial condition u(x,0)=exp(-100*x^2)*cos(pi*x/2) on the bounds
%     -1<x<1,0<t<1 with periodic x boundary conditions. It employs a 3-tier 
%     semi-implicit finite difference scheme, and uses gauss-seidel to solve 
%     the resulting implicit system. Used for problem #4 of the final.
% Inputs:
%     subs:   the spatial step size
%     tol:    approximation tolerance
%     it:     max number of iterations for the gauss-seidel method
%             (interpreted as infinite if 0)
% Outputs:
%     meshx:      meshgrid of x values used for plotting. (created by the meshgrid function)
%     mesht:      meshgrid of t values used for plotting. (created by the meshgrid function)
%     w:          (p+1)x(m+1) matrix of approximations w(i,j) to the solutions at points x_ij
%                 for i=1,..,(p+1), j=1,...,(m+1)
%     frames:     array of frames of the plot of the error at each iteration
%                 of GS, (for error analysis and debugging purposes)
function [meshx,mesht,w,frames]=fdWrap(subs,tol,it)
    %variables (bounds,subintervals,etc.)
    a=-1;b=1;c=0;d=1;dx=subs;dt=dx/2;p=floor((b-a)/dx);m=floor((d-c)/dt);
    e=.1;L=(2*dt)/(dx^2);G=dt/dx;beta=e*L;mu=3+(2*beta);
    x=transpose([a:dx:b]);t=transpose([c:dt:d]);
    
    %we set up a grid matrix with exact values along the t=0 border and 0 
    %everywhere else. 
    [meshx,mesht]=meshgrid(x,t);frames=[];
    w=zeros(p+1,m+1);
    w(:,1)=exp(-100*x.^2).*cos(pi*x/2);
    
    %we fill w(:,2) with our initial approximation 
    mid=round((p+1)/2);
    for i=1:(p+1)
        w(i,2)=g(x(i))-dt*g(x(i))*g1(x(i))+e*dt*g2(x(i))+dt*f(x(i),0);
        %surf(meshx,mesht,w');
    end
    
    %now we execute our finite difference scheme on the environment we set
    %up, using gauss seidel to solve the system
    l=1;norm=[];ni=[];
    set(gcf,'Renderer','zbuffer');
    while true
        for j=2:m
            for k=mid:(mid+p-1)
                modid=mod(k,p);
                i=modid+1;
                if i==1;
                    %approximate w(1,j+1) based on periodic boundaries
                    %(also determines w(p+1,j+1)
                    Cx=-w(i,j-1)+G*w(i,j-1)*(w(i+1,j-1)-w(p,j-1));
                    Bx=4*w(i,j)-2*G*w(i,j)*(w(i+1,j)-w(p,j));
                    ax=beta*(w(i+1,j+1)+w(p,j+1));
                    %calculate next approximation using gauss-seidel
                    z=(1/mu)*(ax+Bx+Cx);
                    %residual vector analysis
                    if norm<abs(z-w(i,j+1))
                        norm=abs(z-w(i,j+1));ni=[i,j+1];
                    end;
                    %store next approximation
                    w(i,j+1)=z;w(p+1,j+1)=z;
                    %surf(meshx,mesht,w');
                elseif i==p
                    %approximate w(p,j+1) based on periodic boundaries
                    Cx=-w(i,j-1)+G*w(i,j-1)*(w(1,j-1)-w(i-1,j-1));
                    Bx=4*w(i,j)-2*G*w(i,j)*(w(1,j)-w(i-1,j));
                    ax=beta*(w(1,j+1)+w(i-1,j+1));
                    %calculate next approximation using gauss-seidel
                    z=(1/mu)*(ax+Bx+Cx);
                    %residual vector analysis
                    if norm<abs(z-w(i,j+1))
                        norm=abs(z-w(i,j+1));ni=[i,j+1];
                    end;
                    %store next approximation
                    w(i,j+1)=z;
                    clearvars i;
                    %surf(meshx,mesht,w'); 
                else
                    %calculate values from boundary conditions and previous time steps
                    Cx=-w(i,j-1)+G*w(i,j-1)*(w(i+1,j-1)-w(i-1,j-1));
                    Bx=4*w(i,j)-2*G*w(i,j)*(w(i+1,j)-w(i-1,j));
                    ax=beta*(w(i+1,j+1)+w(i-1,j+1));
                    %calculate next approximation using gauss-seidel
                    z=(1/mu)*(ax+Bx+Cx);
                    %residual vector analysis
                    if j==2 && i==(mid+1)
                        norm=abs(z-w(i,j+1));ni=[i,j+1];
                    end
                    if norm<abs(z-w(i,j+1))
                        norm=abs(z-w(i,j+1));ni=[i,j+1];
                    end;
                    %store next approximation
                    w(i,j+1)=z;
                    %surf(meshx,mesht,w'); 
                    
                end
            end
        end
        %snippet for error graphing purposes
%         surf(meshx,mesht,w');
%         view(-1.545000000000000e+02,31.999999999999982);
%         axis([-1,1,0,1,-1,1]);
%         drawnow;
%         frames=[frames;getframe];       
        
        %print current norm and current gauss-seidel iteration count
        if l>1
            %disp(' ');
            fprintf(repmat(sprintf('\b'),1,mlen));
        end
        msg=sprintf('l=%d\tnorm=%d\t@\t(%d,%d)',l,norm,ni(1),ni(2));mlen=length(msg+3);
        fprintf(msg);
        
        %check norm against tolerance and iteration count (terminating
        %conditions)
        if norm<=tol || l==it
            disp(' ');
            break;
        end
        
        %increment iteration count
        l=l+1;
    end
    %initial condition u(x,0)=g(x)
    function sol=g(xi)
        sol=exp(-100*xi^2)*cos(pi*xi/2);
    end
	%finite difference approximation to g'(x)
    function sol=g1(xi)
        forward=~boolean(xi-1)*g(x(2))+boolean(xi-1)*g(xi+dx);
        backward=~boolean(xi+1)*g(x(p))+boolean(xi+1)*g(xi-dx);
        sol=(forward-backward)/(2*dx);
    end
    %finite difference approximation to g''(x)
    function sol=g2(xi)
        forward=~boolean(xi-1)*g(x(2))+boolean(xi-1)*g(xi+dx);
        backward=~boolean(xi+1)*g(x(p))+boolean(xi+1)*g(xi-dx);
        sol=(forward-2*g(xi)+backward)/(dx^2);
    end
end