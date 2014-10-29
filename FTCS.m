% Joey Willhite
% Description:
%     Function to approximate solutions to the heat equation in two dimensions
%     with u(0,0,t)=u(using FTCS difference approximation
% Inputs:
%     f:            initial condition u(x,y,0)=f(x,y)
%     intervals:    array [X,Y,T] for intervals 0<x<X, 0<y<Y, 0<t<T
%     bounds:       array [g1,g2,h1,h2] of boundary conditions u(0,y,t)=g1, u(X,y,t)=g2
%                       u(x,0,t)=h1, u(x,Y,t)=h2
%     coeffs:       array [a,b] of coefficient a and b for the PDE u_tt=au_xx+bu_yy
%     params:       array [dx,dy,dt] of x,y,t parameters respectively
%     **NOTE**
%         for dx,dy,dt: if <=1 will be treated as step size, if >1 will be treated
%         as number of subintervals
% Outputs:
%     w:      m*n*r matrix w approximating u(x_i, y_j, t_s) for i=1..m, j=1...n, s=1...r
function w=FTCS(f,intervals,bounds,coeffs,params)
    %pull our boundary conditions and coefficients from our inputs
    bounds=num2cell(bounds); [g1,g2,h1,h2]=deal(bounds{:});
    coeffs=num2cell(coeffs); [a,b]=deal(coeffs{:});
    g1=symfun(g1,symvar(g1)); g2=symfun(g2,symvar(g2));
    h2=symfun(h2,symvar(h2)); h1=symfun(h1,symvar(h1));
    
    %loop to set our proper x,y,t intervals and indices
    subint=[0,0,0];
    for i=1:length(params)
        %if the inputed parameter is less than 1
        if params(i)<1
            %it will be treated as the step size, and m will be calculated
            subint(i)=floor(intervals(i)/params(i));
        else
            %otherwise it will be treated as the subinterval count, in which
            %case the step size and subinterval count are adjusted accordingly
            subint(i)=params(i); params(i)=intervals(i)/subint(i);
        end
    end
    intervals=num2cell(intervals); [X,Y,T]=deal(intervals{:});
    params=num2cell(params); [dx,dy,dt]=deal(params{:});
    subint=num2cell(subint); [m,n,r]=deal(subint{:});
    alpha=(a*dt)/(dx^2); beta=(b*dt)/(dy^2);
    
    %create our matrix w, and fill in w(x,y,0) with initial conditions from f 
    w=zeros(m+1,n+1,r+1);
    for i=1:m+1
        for j=1:n+1
            w(i,j,1)=f((i-1)*dx,(j-1)*dy);
        end
    end
    %fill w with boundary values from g1,g2,h1,h2
    for i=1:r
        %boundary u(0,y,t)=g1
        w(1,:,(i+1))=h1(i*dt,linspace(0,X,m+1));
        %boundary u(X,y,t)=g2
        w(m+1,:,(i+1))=h2(i*dt,linspace(0,X,m+1));
        
        %boundary u(x,0,t)=h1
        w(:,1,(i+1))=g1(i*dt,linspace(0,Y,n+1));
        %boundary u(x,Y,t)=h2
        w(:,n+1,(i+1))=g2(i*dt,linspace(0,Y,n+1));
        
    end
    
    %loop to fill w with approximations
    for i=1:r
        for j=2:m
            for k=2:n
            w(j,k,i+1)=(1-2*alpha-2*beta)*w(j,k,i)+beta*(w(j-1,k,i)+w(j+1,k,i))+alpha*(w(j,k-1,i)+w(j,k+1,i));
            end
        end
    end
        
end