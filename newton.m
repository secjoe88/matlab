function f = newton(symExpr, pNot, tol)

%turn the inputed symbolic expression into an inline
%function and create an inline function for its derivative
func = inline(char(symExpr));
funcPrime = inline(diff(symExpr));

%if the initial approximation is a root, return it
if func(pNot)==0
	f=pNot;

%otherwise find next closest approximation
else
	p=pNot-(func(pNot)/funcPrime(pNot));
end

%if the next closest approximation is within tol of the current one 
%(which includes it being a root itself), return it
if (abs(p-pNot)<=tol) 
	f=p;

%otherwise recurse
else
	newton(symExpr, p, tol) 
end
