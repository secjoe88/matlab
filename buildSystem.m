% Description:
%     a script to build the system Ax=b for the elliptic pde u_xx+u_yy=0 with boundary
%     condition g(x,y)=log(x^2+y^2), using the implicit FTCS scheme as an approximation
%     to the PDE, at with inputed dx.
% inputs:
%     delta:      timestep interval dx=dy     
% outputs:
%     A:          (n-1)(m-1)x(n-1)(m-1) matrix A for Ax=b
%     b:          vector b for Ax=b
function [A,b]=buildSystem(delta)
    m=floor(1/delta);
    %construct blocks, and form system as block matrix
    dblock=4*eye(m-1)-diag(ones(1,m-2),1)-diag(ones(1,m-2),-1);
    oblock=-eye(m-1);
    A=kron(eye(m-1),dblock)+kron(diag(ones(1,m-2),1),oblock)+kron(diag(ones(1,m-2),-1),oblock);
    disp(['Constructed A for dx=dy=',num2str(delta)]);
    %calculate b based on boundary conditions
    xs=1+delta:delta:2-delta; ys=delta:delta:1-delta;
    b=zeros(1,(m-1)^2);
    for i=1:m-1
        for j=1:m-1
            p=i+(m-1-j)*(m-1);
            if i==1
                b(p)=b(p)+g(1,ys(j));
            end
            if i==(m-1)
                b(p)=b(p)+g(2,ys(j));
            end
            if j==1
                b(p)=b(p)+g(xs(i),0);
            end
            if j==(m-1)
                b(p)=b(p)+g(xs(i),1);
            end
        end
    end
    disp(['Constructed b for dx=dy=',num2str(delta)]);
    function sol=g(x,y)
        sol=log(x^2+y^2);
    end
end
