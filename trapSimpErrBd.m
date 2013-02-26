function [trapErrBd, simpErrBd]=trapSimpErrBd(a,b,f)
    h=b-a;
    syms('f2(x)');f2(x)=diff(f,2);
    syms('f4(x)');f4(x)=diff(f,4);
    
    trapErrBd=-((h^3)/12)*(f2(a)-f2(b));
    simpErrBd=-((h^5)/90)*(f4(a)-f4(b));
    
end