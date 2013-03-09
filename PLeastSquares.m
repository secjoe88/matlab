% Joey Willhite
% Description:
%     A funciton to create polynomial approximations of degree n to arbitrary
%     datasets or to a function using the Least Squares approximation method.
% Inputs:
%     option: String representing whether the method is to approximate a collection
%             of data (with option='data') or a function (with
%             option='function')
%     n: The degree of the desired polynomial approximation
%     *Note*
%       Function will ask for user input for data or function/interval
%       depending on value of 'option' parameter.
% Outputs:
%     solutions: Struct containing the calculated coefficients of the
%     polynomial approximation
function solutions=PLeastSquares(option, n)
    a_i=sym('a', [1 n+1]);
    eqs=cell(n+1, 1);
    if strcmp(option,'data')==1 
        x_i=input('Input x values: ');
        y_i=input('Input y values: ');
        for i=1:n+1
            eqs{i}=buildEq(option, x_i, y_i, a_i, n, i);
        end
    elseif strcmp(option, 'function')==1
        syms x;
        func=input('Input function (as symbolic function of x): ');
        interval=zeros(1,2);
        interval(1)=input('Beginning of interval: ');
        interval(2)=input('End of interval: ');
        for i=1:n+1
            eqs{i}=buildEq(option, func, interval, a_i, n, i);
        end
        
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