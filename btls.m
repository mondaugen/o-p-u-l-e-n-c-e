function [s]=btls(x,dx,f,df,al,be,fopts)
    % Backtracking line search
    s=1;
    while f(x+s*dx,fopts)>(f(x,fopts)+al*s*df(x,fopts)'*dx)
        s=be*s;
    end
end
