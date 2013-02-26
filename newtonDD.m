% Joey Willhite
% Inputs:
%     x=array of x values {x_0,...,x_n}
%     F=array of f(x_n) values {f(x_0),...,f(x_n)}
% Outputs:
%     coefficients=array of coefficients for the terms of the 
%        interpolating polynomial
% *Note*
%     Before the program returns the coefficents of the reulting 
%     polynmial, it prints out the entire divided difference chart


function [coefficients]= newtonDD(x, F)
    fLength=size(F);
    fLength=fLength(2);
    chart=zeros(fLength);
    chart(:,1)=F;
    coefficients=zeros(1, 4);
    coefficients(1,1)=F(1,1);
    for i=2:fLength,
        for j=2:i,
            chart(i,j)=(chart(i,j-1)-chart(i-1, j-1))/(x(i)-x(i-(j-1)));
            if i==j,
                coefficients(1,i) = chart(i,j);
            end
        end
    end
    disp(chart);
end


 