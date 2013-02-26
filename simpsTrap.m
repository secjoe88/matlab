% Joey Willhite
% Description:
%     a simple function to use Simpson's Rule and the Trapezoid Rule
%     to evaluate the integral of a function over a given interval
% Inputs:
%     a: beginning of interval
%     b: end of interval
%     f: function (as symexpr)
% Outputs:
%     simRes: function evaluated with Simpson's Rule
%     trapRes: function evaluated with Trapezoid Rule

function [simRes, trapRes]=simpsTrap(a,b,f)
    h=b-a;
    trapRes=(h/2)*(f(a)+f(b));
    simRes=(h/3)*(f(a)+f((a+b)/2)+f(b));
end