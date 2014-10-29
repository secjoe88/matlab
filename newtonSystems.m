% Joey Willhite
% Description:
%     a function to approximate the solution to a nonlinear system of equations
%     using newton's method.
% Inputs: 
%     functions:  a (column)cell array of functions to solve (as symfunc)
%     init:       an (row)array of initial values
%     n:          number of functions/unknowns to solve for
%     terminator: terminating case. if terminator>1 it will be treated as an 
%                 iteration. if terminator<1 it will be treated as a tolerance
% Output:
%     solutions:  Sequence of approximations converging to solution

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
        solutions=[next; newtonSystems(functions, next, n, terminator)];
    else
        %if terminating case is iterative, we recurse
        solutions=[next;newtonSystems(functions, next, n, terminator-1)];
    end

    
