% Joey Willhite
% Description:
%     A function to use the Poisson equation finite difference algorithm (12.4 Burden & Faires)
%     optimized with SOR to approximate solutions to elliptic PDEs. 
% Inputs:
%     **NOTE**
%     Helper functions f(x,y) and g(x,y) for the functions u_xx+u_yy=f(x,y) and boundary
%     condition g(x,y) must be coded explicitly in the helper function below. This is
%     necessary to allow for scaling to large systems
%     xbound:     Array [a,b] of x bounds
%     ybound:     Array [c,d] of y bounds
%     nh:         If nh>1 it will be treated as number of subintervals, if
%                   nh<=1 it will be treated as dx
%     mk:         If mk>1 it will be treated as number of subintervals, if
%                   mk<=1 it will be treated as dy
%     tol:        Tolerance
%     N:          Max iterations
%                 **NOTE** if N=0, will be treated as infinite
% Outputs:
%     x:          Array [x_1, x_2,...,x_n-1] for x_i=a+ih
%     y:          Array [y_1,y_2,...,y_m-1] for y_j=c+jk
%     w:          (n-1)*(m-1) matrix of approximations w_ij to interior points u(x_j,y_i)
%     frames:     array of graph frames showing evolution of error over time
function [x,y,w,frames]=poissonFD(xbound, ybound, nh, mk, tol, N)
    %steps 1 adjusted for robustness
    if(nh>1)
        n=nh;
        h=(xbound(2)-xbound(1))/n;
    else
        h=nh;
        n=floor((xbound(2)-xbound(1))/h);
    end
    if(mk>1)
        m=mk;
        k=(ybound(2)-ybound(1))/m;
    else
        k=mk;
        m=floor((ybound(2)-ybound(1))/k);
    end
    %define omega for SOR method
    omega=4/(2+sqrt(4-(cos(pi/m)+cos(pi/n))^2));
    %steps 2-5
    x=0; y=0; w=zeros(n-1,m-1); a=xbound(1);b=xbound(2); c=ybound(1);d=ybound(2); 
    for i=1:(n-1)
        x(i)=a+i*h;
    end
    for i=1:(m-1)
        y(i)=c+i*k;
    end
    lambda=h^2/k^2; mu=2*(1+lambda);l=1;
    %step 6
    frames=[];mlen=0;
    while true
        %step 7 (7-9: Calculate top row and top corners)
        z=(1-omega)*w(1,m-1)+omega*(-h^2*f(x(1), y(m-1))+g(a, y(m-1))+lambda*g(x(1),d)+lambda*w(1,m-2)+w(2,m-1))/mu;
        norm=abs(z-w(1,m-1));
        w(1,m-1)=z;
        %step 8
        for i=2:(n-2)
            z=(1-omega)*w(i,m-1)+omega*(-h^2*f(x(i),y(m-1))+lambda*g(x(i),d)+w(i-1,m-1)+w(i+1,m-1)+lambda*w(i,m-2))/mu;
            if abs(w(i,m-1)-z)>norm
                norm=abs(w(i,m-1)-z);
            end
            w(i,m-1)=z;
        end
        %step9
        z=(1-omega)*w(n-1,m-1)+omega*(-h^2*f(x(n-1),y(m-1))+g(b,y(m-1))+lambda*g(x(n-1),d)+w(n-2,m-1)+lambda*w(n-1,m-2))/mu;
        if abs(w(n-1,m-1)-z)>norm
            norm=abs(w(n-1,m-1)-z);
        end
        w(n-1,m-1)=z;
        %step 10 (10-13: Calculate middle rows and edges)
        for j=(m-2):-1:2
            %step 11
            z=(1-omega)*w(1,j)+omega*(-h^2*f(x(1),y(j))+g(a,y(j))+lambda*w(1,j+1)+lambda*w(1,j-1)+w(2,j))/mu;
            if abs(w(1,j)-z)>norm
                norm=abs(w(1,j)-z);
            end
            w(1,j)=z;
            %step 12
            for i=2:(n-2)
                z=(1-omega)*w(i,j)+omega*(-h^2*f(x(i),y(j))+w(i-1,j)+lambda*w(i,j+1)+w(i+1,j)+lambda*w(i,j-1))/mu;
                if abs(w(i,j)-z)>norm
                    norm=abs(w(i,j)-z);
                end
                w(i,j)=z;
            end
            %step 13
            z=(1-omega)*w(n-1,j)+omega*(-h^2*f(x(n-1),y(j))+g(b,y(j))+w(n-2,j)+lambda*w(n-1,j+1)+lambda*w(n-1,j-1))/mu;
            if abs(w(n-1,j)-z)>norm
                norm=abs(w(n-1,j)-z);
            end
            w(n-1,j)=z;
        end
        %step 14 (14-16: Calculate bottom row and bottom corners)
        z=(1-omega)*w(1,1)+omega*(-h^2*f(x(1),y(1))+g(a,y(1))+lambda*g(x(1),c)+lambda*w(1,2)+w(2,1))/mu;
        if abs(w(1,1)-z)>norm
            norm=abs(w(1,1)-z);
        end
        w(1,1)=z;
        %step 15
        for i=2:(n-2)
            z=(1-omega)*w(i,1)+omega*(-h^2*f(x(i),y(1))+lambda*g(x(i),c)+w(i-1,1)+lambda*w(i,2)+w(i+1,1))/mu;
            if abs(w(i,1)-z)>norm
                norm=abs(w(i,1)-z);
            end
            w(i,1)=z;
        end
        %step 16
        z=(1-omega)*w(n-1,1)+omega*(-h^2*f(x(n-1),y(1))+g(b,y(1))+lambda*g(x(n-1),c)+w(n-2,1)+lambda*w(n-1,2))/mu;
        if abs(w(n-1,1)-z)>norm
            norm=abs(w(n-1,1)-z);
        end
        w(n-1,1)=z;
        %added bit for error graphing purposes
%         [xs,ys]=meshgrid(x,y);
%         exact=log(xs.^2+ys.^2);
%         surf(xs,ys,abs(transpose(exact)-w));
%         axis([1,2,0,1,0,1.5]);
%         set(gcf,'Renderer', 'zbuffer');
%         frames=[frames;getframe()];
        if l>1
            fprintf(repmat(sprintf('\b'),1,mlen));
        end
        msg=sprintf('l=%d\tnorm=%d',l,norm);mlen=length(msg+1);
        fprintf(msg);
        %step 17-19
        if norm<=tol || l==N
            break;
        end
        %step 20
        l=l+1;
    end
    disp(' ');
    %step 21
    if l==n
        display('Max iterations Exceeded');
    end
    
    function sol=f(x,y)
        sol=0;
    end
    function sol=g(x,y)
        sol=log(x^2+y^2);
    end
end
