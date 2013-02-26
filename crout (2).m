% Joey Willhite
% Description:
%     A function to use Crout's Algorithm to solve tridiagonal linear
%     systems.
% Inputs:
%     a: A vector containing the entries of the middle diagonal of the 
%        coefficient matrix of the linear system to be solved.
%     b: A vector containing the entries of the bottom diagonal of the 
%        coefficient matrix of the linear system to be solved.
%     c: A vector containing the entries of the top diagonal of the 
%        coefficient matrix of the linear system to be solved.
%     s: A vector containing the (n+1)th column of the coefficient matrix
%        of the linear system to be solved.
%     n: The number of unknowns in the linear system to be solved for
% Outputs:
%     This function has no outputs.
%     *Note* the function prints out the solutions to the linear system
function Crout(a, b, c, s, n)
    %instantiate l, u, z, x vectors
    l=zeros(1, n);
    u=zeros(1, n-1);
    z=zeros(1, n);
    x=zeros(1, n);
    %Step 1
    l(1)=a(1);
    u(1)=(c(1)/l(1));
    z(1)=s(1)/l(1);
    %Step 2
    for i=2:(n-1)
        l(i)=a(i)-(b(i-1)*u(i-1));
        u(i)=(c(i)/l(i));
        z(i)=(s(i)-(b(i-1)*z(i-1)))/l(i);
    end
    %Step 3
    l(n)=a(n)-(b(n-1)*u(n-1));
    z(n)=(s(n)-(b(n-1)*z(n-1)))/l(n);
    %Step 4
    x(n)=z(n);
    %Step 5
    for i=(n-1):-1:1
        x(i)=z(i)-(u(i)*x(i+1));
    end
    
    %Step 6 (also displays L and U matrices
    disp('L=');
    disp(num2str(diag(l)+diag(b,-1)));
    disp('U=');
    disp(num2str(eye(n)+diag(u, 1)));
    for i=1:(n-1)
        fprintf(['x_' num2str(i) '=' num2str(x(1)) ', ']);
    end
    disp(['x_' num2str(n) '=' num2str(x(n))]);
end