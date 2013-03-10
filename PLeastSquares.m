% Joey Willhite
% Description:
%     A funciton to create polynomial approximations of degree n to arbitrary
%     datasets or to a function using the Least Squares approximation method.
% Inputs:
%     option: String representing whether the method is to approximate a collection
%             of data (with option='data') or a function (with
%             option='function')
%     arg1: PLeastSquares expects a function if option value is 'function'
%       or x-values of data if option value is 'data'
%     arg2: PLeastSquare expects an array representing the interval of
%       approximation if option value is 'function' or y-values of data if
%       option value is 'data'
%     n: The degree of the desired polynomial approximation
%Outputs:
%     solutions: Struct containing the calculated coefficients of the
%     polynomial approximation
function solutions=PLeastSquares(option, arg1, arg2, n)
    a_i=sym('a', [1 n+1]);
    eqs=cell(n+1, 1);
    for i=1:n+1
            eqs{i}=buildEq(option, arg1, arg2, a_i, n, i);
    end
    solutions=solve([eqs{1:length(eqs)}]);
end

%Function to build normal equations:
%   This function builds the corresponding normal equations to be solved
%   for the coefficients of the polynomial approximation. It uses a switch
%   statement based on the 'option' parameter to determine whether it will
%   build normal equations pursuant to an approximation of data, with arg1
%   and arg2 being the x and y values respectively, or pursuant to an
%   approximation of a function, with arg1 and arg2 being the function and
%   interval of aproximation respectively.
function eq=buildEq(option, arg1, arg2, coeffs, degree, eq_no)
    eq=0;
    j=1;
    switch option
        case 'data'
            for i=(eq_no-1):degree+(eq_no-1)
                eq=eq+(coeffs(j)*sigmaSingle(arg1,i));
                j=j+1;
            end
            eq=eq==sigmaDouble(arg1,arg2,eq_no-1);
        case 'function'
            syms x;
            for i=(eq_no-1):degree+(eq_no-1)
                eq=eq+(coeffs(j)*int(x^i, arg2(1), arg2(2)));
                j=j+1;
            end
            eq=eq==int((x^(eq_no-1))*arg1,arg2(1),arg2(2));
    end
    
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