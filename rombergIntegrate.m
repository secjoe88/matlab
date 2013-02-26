% Joey Willhite
% Description:
%     A function to use the Romberg technique to approximate
%     an integral over a finite region-also can define desired 
%     precision and max num of intervals
% Inputs:
%     a: beginning of interval(x_0)
%     b: end of interval(x_n)
%     f: the integrand
%     n: max number of subintervals (2^(n-1))
%     tol: maximum precision
% Output:
%     Prints out the full Romberg tree one line at a time followed 
%     by the estimated precision of the final approximation

function rombergIntegrate(a, b, f, n, tol)
    rombergTree=zeros(2,n);
    %Step 1
    h=b-a;
    rombergTree(1,1)=(h/2)*(f(a)+f(b));
    %Step 2
    disp(rombergTree(1,1));
    %Step 3
    for i=2:n,
        %Step 4
        rombergTree(2,1)=(1/2)*(rombergTree(1,1)+trapezoid(i,h));
        %Step 5
        for j=2:i,
            rombergTree(2,j)=rombergTree(2,j-1)+...
                ((rombergTree(2,j-1)-rombergTree(1,j-1))/...
                (4^(j-1)-1));
        end
        %Step 6
        disp(rombergTree(2,:));
        %Tolerance check-if within tolerance then print
        %tolerance and break
        if abs(rombergTree(1,i-1)-rombergTree(2,i))<=tol
            format long; 
            str=[num2str(rombergTree(1,i-1), 10), ...
                 '-', num2str(rombergTree(2,i), 10), '=', ...
                num2str(rombergTree(1,i-1)-rombergTree(2,i), 10)];
            disp(str);
            format short;
            break
        else
            %Otherwise proceed to Step 7
            h=h/2;
            %Step 8
            rombergTree(1,:)=rombergTree(2,:);
        end
    end
    
    %Function to evaluate the Trapezoid rule implemented in
    %Step 4
    function value=trapezoid(i,h)
        value=0;
        for k=1:2^(i-2),
            value=value+h*(f(a+(k-0.5)*h));
        end
    end
end

