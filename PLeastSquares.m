% Joey Willhite
% Description:
%     A funciton to create polynomial approximations of degree n to arbitrary
%     datasets using the Least Squares approximation method.
% Inputs:
%     x_i: The x values of the given dataset
%     y_i: The y values of the given dataset
%     a_i: Symbolic variables corresponding to the coefficients of the
%          approximating polynomial
%     n: The degree of the desired polynomial approximation
% Outputs:
%     solutions: Struct containing the calculated coefficients of the
%     polynomial approximation
function solutions=PLeastSquares(x_i, y_i, a_i, n)
    eqs=cell(n+1, 1);
    for i=1:n+1
        eqs{i}=buildEq(x_i,y_i,a_i,n,i);
    end
    solutions=solve([eqs{1:length(eqs)}]);
end

%Build normal equations
function eq=buildEq(xs, ys, coeffs, degree, eq_no)
    eq=0;
    j=1;
    for i=(eq_no-1):degree+(eq_no-1)
        eq=eq+(coeffs(j)*sigmaSingle(xs,i));
        j=j+1;
    end
    eq=eq==sigmaDouble(xs,ys,eq_no-1);
end

%compute sigma function for the sum of x values to the degree specified
function sigma=sigmaSingle(xs,degree)
    sigma=0;
    for i=1:length(xs)
        sigma=sigma+(xs(i)^degree);
    end
end

%compute sigma function for sum of y*x^degree
function sigma=sigmaDouble(xs, ys, degree)
    sigma=0;
    for i=1:length(xs)
        sigma=sigma+(ys(i)*(xs(i)^degree));
    end
end