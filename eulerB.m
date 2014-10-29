% Description:
%     a function to solve systems of first order differential equations using
%     the backward Euler's method(resulting nonlinear sytems solved using
%     Newton's method).
% Inputs:
%     functions:      cell (column) array of functions fn(t,y1,..,yn) for the system
%                     (y1', y2',...,yn')=(f1, f2,...,fn)
%     domain:         domain [a,b] of approximation 
%     init:           array [y1(a), y2(a),...,yn(a)] of initial values
%     h:              step size
% Output:
%     approximation:  (column) array of approximations [tn, y(tn), y'(tn)]
function approximation=eulerB(functions, domain, init, h)
    %known w values
    ws=transpose(init);
    %variables we'll be using to solve the system later
    w_next=sym('w',[length(functions),1]);
    
    %we set up our functions to solve as a linear system
    temp_funcs={};
    
    %begin algorithm loop
    a=domain(1);
    while a<domain(2)
        for i=1:length(functions)
            temp_funcs{i}=vpa(w_next(i)-(ws(i,length(ws(i,:)))+h*subs(functions(i),symvar([functions{:}],15),transpose(w_next))));
        end
        %solve linear system using newton's method (find next iteration)
        res=newtonSystems(temp_funcs,zeros(1,length(temp_funcs)),length(temp_funcs), 10^-5);
        ws=[ws,transpose(res)];
        
        
       %increment count
        a=a+h;
    end
    approximation=[transpose(domain(1):h:domain(2)),transpose(ws)];
end

function solutions = newtonSystems(functions, init, n, terminator)
    %%terminating condition (for iterative terminator)
    if terminator==0
        solutions=init
        return
    end
    %evaluate our functions at our initial guesses init
    Fx=zeros(n,1);
    for i=1:n
        Fx(i)=subs(functions{i},symvar([functions{:}]), init);
    end
    
    %%check if initial guess is indeed a solution.
    if Fx==0
        solutions=init
        return
    end
    
    %calculate the jacobian, and evaluate it at our initial guess
    Jacobian=jacobian(functions, symvar([functions{:}]));

    %and evaluate it at our initial guess as well
    Jacobian=double(subs(Jacobian, symvar([functions{:}]), init));

    %solve the linear system J(x)y=-F(x), and create our next approximatoin
    y=double(Jacobian\(-Fx));
    next=init+transpose(y);
    
    
    %if the terminating mechanism is a tolerance
    if terminator<1
        %return if approximation is within tolerance
        if norm(transpose(y))<terminator
            solutions=next;
            return
        end
        %otherwise recurse with same tolerance
        solutions=newtonSystems(functions, next, n, terminator);
    else
        %if terminating case is iterative, we recurse
        solutions=newtonSystems(functions, next, n, terminator-1);
    end
end