function gaussianElimination(A, n)
    A = eliminate(A, n);
    
    
    
end

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
function x =backwardSub(A, n)
    x=zeros(1, n);
    x(1,n)=A(n, n+1)/A(n,n);
end
            
        