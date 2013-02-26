% Joey Willhite
% Description:
%     simple function to evaluate the integral of a function on a given
%     interval using the Trapezoid, Simpson's, and Midpoint Composite 
%     Rules
% Inputs:
%     a: beginning of interval
%     b: end of interval
%     n: number of subintervals
%     f: function to be evaluated
% Outputs:
%     xT: result obtained with Trapezoid Rule
%     xS: result obained with Simpson's Ruls
%     xM: result obtained with Midpoint Rule


function [xT, xS, xM] = compTSMRule(a,b,n,f)
    %Step 1
    x=0;
    h=(b-a)/n;
    
    %Step 2
    xi0=f(a)+f(b);
    xi1=0;
    xi2=0;
    
    %Step 3-5
    for i=1:n-1,
        x=a+i*h;
        if (mod(i,2)==0)
            xi2=xi2+f(x);
        else
            xi1=xi1+f(x);
        end
    end
    
    %Step 6
    xS=h*(xi0+(2*xi2)+(4*xi1))/3;
    xT=h*(xi0+2*(xi2+xi1))/2;
    xM=2*h*(xi1);
        