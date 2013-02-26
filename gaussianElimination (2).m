% Joey Willhite
% Description:
%     Function to calculate the solutions to a system of 
%     linear equations using gaussian elimination with backward
%     substitution
% Inputs:
%     A: nxn+1 Augmented matrix representing the linear system
%     n: number of unknowns for the system
% Ouputs:
%     The function has no outputs
%     *Note* function prints out solutions to unknowns line by 
%     line

function gaussianElimination(A, n)
    A = eliminate(A, n);
    x=backwardSub(A, n);
    for i=1:n
        str=['x_', num2str(i), '=', num2str(x(1, i))];
        disp(str);
    end
    
    
    
end
%function to perform elimination step
function A = eliminate(A, n)
    p=0;
    %Step 1: Begin elimination process
    for i=1:(n-1)
        %Step 2: find smallest p such that a_pi!=0
        for p=i:n
            if A(p,i)~=0
                break;
            end
            %if there is no p, output that the shit is sketch
            if p==n
                error('No unique solution');
            end
        end
        %Step 3: if a_pi is not on the diagonal, swap rows i and p
        if p~=i
            tempRow = A(p,:);
            A(p,:)=A(i,:);
            A(i,:)=tempRow;
        end
        %Step 4: row reduction step
        for j=(i+1):n
            %Step 5:set fraction value
            m_ji=A(j,i)/A(i,i);
            %Step 6: Reduce row
            A(j,:)=(A(j,:)-(m_ji*A(i,:)));
        end
    end
    %Step 7: Check if elimination step reduced to system with non-unique
        %solution
    if A(n, n)==0
        error('No unique solution');
    end
%end elimination step    
end
%function to perform backwards substitution step
function x =backwardSub(A, n)
    x=zeros(1, n);
    %Step 8: Begin backwards substitution
    x(1,n)=A(n, n+1)/A(n,n);
    %Step 9: calculate x_i values
    for i=(n-1):-1:1
        sigma=0;
        %calculate sigma(a_ij*x_j)
        for j=(i+1):n
            sigma=sigma+(A(i,j)*x(1,j));
        end
        x(1,i)=(A(i, n+1)-sigma)/A(i,i);
    end
end
            
        