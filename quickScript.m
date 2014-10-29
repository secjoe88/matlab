% Description:
%     to discern the convergence of the finite difference scheme by approximating
%     the solution to the linear system Ax=b using gaussian elimination for several
%     choices of dx=dy
% inputs:
%     deltas:     array [d1,d2,...,dn] of delta values to test
% outputs:
%     errors:     errors [e1,e2,...,en] where e=max(u_ij-u(x_i,y_j)) for each delta
function errors=quickScript(deltas,g)
    
    k=1;errors=[];
    while k<=length(deltas)
        m=floor(1/deltas(k));
        %construct blocks, and form system as block matrix
        dblock=4*eye(m-1)-diag(ones(1,m-2),1)-diag(ones(1,m-2),-1);
        oblock=-eye(m-1);
        A=kron(eye(m-1),dblock)+kron(diag(ones(1,m-2),1),oblock)+kron(diag(ones(1,m-2),-1),oblock);
        disp(['Constructed A for dx=dy=',num2str(deltas(k))]);
        %calculate b based on boundary conditions
        xs=1+deltas(k):deltas(k):2-deltas(k); ys=deltas(k):deltas(k):1-deltas(k);
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
        disp(['Constructed b for dx=dy=',num2str(deltas(k))]);
        %solve system using gaussian elimination, and store approximation
        w=A\transpose(b);approx=[];
        for i=1:(m-1)
            approx(m-i,:)=w((m-1)*(i-1)+1:i*(m-1));
        end
        disp(['Calculated approximation for dx=dy=',num2str(deltas(k))]);
        
        %create meshgrid of x,y pairs and calculate exact value of
        %log(x^2+y^2)
        [xs,ys]=meshgrid(xs,ys);
        exact=log(xs.^2+ys.^2);
        disp(['Calculated exact solutions for dx=dy=',num2str(deltas(k))]);
        %and store error between approximation and actual solution
        errors=[errors; max(max(abs(approx-exact)))];
        disp(['Error stored for dx=dy=',num2str(deltas(k))]);
        k=k+1;
    end
end
