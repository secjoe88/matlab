function k1 =test(x,y,f)
    h=(x-y)/2;
    a=x;
    b=y;
    disp(int2str(a));
    k1=eval(h*f(a,b));
    disp(int2str(k1));
    k2=k1+1;
end