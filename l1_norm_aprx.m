function [z,p]=l1_norm_aprx(A,b,maxiter=100)
    % L1_NORM_APRX
    %
    % solve the optimization problem
    % min || A*x - b ||_1
    % where ||...||1 is the L1 norm
    %
    % To do this we solve the linear program
    % min sum(y)
    % s.t. -y <= A*x - b <= y
    %
    % We can use the barrier method to solve this. The cost function with barrier
    % parameterized by t is
    % C(x,y) = t*sum(y) - sum(log(y - A*x + b) + log(A*x + y - b))
    %
    % Returns z=[x;y];
    DEBUG=1;

    [M,N]=size(A);
    % choose arbitrary x, make feasible y;
    x=zeros(N,1);
    y=(max(abs(b))+1)*ones(M,1);
    z=[x;y];
    % Barrier parameter increment coefficient
    mu=15;
    %mu=50;
    % We assume that a good approximation can be made, that is A*x_ ~ b and so
    % we make t0 = M/||(A*x)||_1 so that initial M/t approximately same order as A*x_
    %t=M/(norm(A*x,1));
    t=mu*M/norm(y,1);
    %t=100;
    ep=1e-6;
    % Line search parameters
    al=0.15;
    be=0.8;

    for it=1:maxiter
        z=l1_norm_aprx_ncs(z,ep,al,be,t,maxiter);
        x=z(1:N);
        y=z(N+1:N+M);
        if (DEBUG)
            display(sprintf('outer iteration: %d',it));
            display(sprintf('Current cost: %f',sum(y)));
        end
        if (4*M)/t < ep
            break
        end
        t=mu*t;
    end

    fopts{1}=t;
    p=f(z,fopts);

    function [p]=f(z,fopts)
        % Evaluate cost function
        t=fopts{1};
        x=z(1:N);
        y=z(N+1:N+M);
        p=t*sum(y)-sum(log(y-A*x+b)+log(A*x+y-b));
    end

    %function [dp]=df(z,fopts)
    %    % Evaluate gradient of cost function
    %    t=fopts{1};
    %    x=z(1:N);
    %    y=z(N+1:N+M);
    %    tmp1=(y-A*x+b).^(-1);
    %    tmp2=(A*x+y-b).^(-1);
    %    g1=tmp1-tmp2;
    %    g2=t*ones(M,1)-tmp1-tmp2;
    %    dp=[A'*g1;g2];
    %end

    function [s]=l1_norm_aprx_btls(z,dz,dp,al,be,fopts,maxiter=1000)
        % Backtracking line search
        % x is current x
        % dx is current delta x (computed from newton step)
        % dp is current gradient
        % al and be are line search params
        % fopts are passed to the cost function when evaluating

        % put x in a reasonable domain?
        s=1;
        for it=1:maxiter
            z_=z+s*dz;
            x_=z_(1:N);
            y_=z_(N+1:N+M);
            if all((y_-A*x_+b)>=0) && all((A*x_+y_-b)>=0)
                break;
            end
            s*=be;
        end
        p=f(z,fopts);
        while (f(z+s*dz,fopts)>(p+al*s*dp'*dz)) && (s > 0.01)
            s=be*s;
        end
        if (s == 0)
            error('Step size is 0.');
        end
        %s=0.1;
        %q=-al*dp'*dx;
        %p=f(x,fopts);
        %while (p - f(x+s*dx,fopts)) < s*q
        %    s=s*be;
        %end
    end

    function [z]=l1_norm_aprx_ncs(z,ep,al,be,t,maxiter=100)
        % Newton centering step
        fopts{1}=t;
        for it=1:maxiter
            x=z(1:N);
            y=z(N+1:N+M);
            tmp1=y-A*x+b;
            tmp2=A*x+y-b;
            D1=diag(tmp1)^(-2);
            D2=diag(tmp2)^(-2);
            tmp1=tmp1.^(-1);
            tmp2=tmp2.^(-1);
            g1=tmp1-tmp2;
            g2=t*ones(M,1)-tmp1-tmp2;
            %D=4*D1*D2*(D1+D2)^(-1);
            D=2*(diag(y)^2+diag(b-A*x)^2)^(-1);
            g=g1+(D1-D2)*(D1+D2)^(-1)*g2;
            % improve stability of solve
            D_=D^(1/2);
            A_=D_*A;
            g_=D^(-1/2)*g;
%            dx=linsolve(A'*D*A,-A'*g);
            dx=ols(-g_,A_);
            dy=(D1+D2)^(-1)*((D1-D2)*A*dx-g2);
            dz=[dx;dy];
            % evaluate gradient
            dp=[A'*g1;g2];
            lamb2=-dp'*dz;
            if (DEBUG)
                display(sprintf('inner iteration: %d',it));
                display(sprintf('lambda^2: %f',lamb2));
            end
            if (lamb2/2<=ep)
                break;
            end
            s=l1_norm_aprx_btls(z,dz,dp,al,be,fopts);
            %s=0.05;
            if (DEBUG)
                display(sprintf('line search result: %f',s));
            end
            z=z+s*dz;
        end
    end
end
