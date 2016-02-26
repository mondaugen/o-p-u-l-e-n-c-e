function [z]=bpdn(A,b,la,maxiter=100)
    % BPDN
    %
    % Basis pursuit denoising.
    % Solve the optimization problem
    % min ||A*x - b||_2 + la*u
    % s.t. -u < x < u
    %
    % la is the parameter that controls sparsity of the solution vector. la near 0
    % does not encourage sparsity. Greater la does.
    %
    % Returns z=[x;u] where x is the solution and u is the variable boundary of x.
    DEBUG=0;

    [M,N]=size(A);
    % choose arbitrary x, make feasible u
    x=zeros(N,1);
    u=2*ones(N,1);
    z=[x;u;];
    % Barrier parameter increment coefficient
    mu=1.1;
    % We assume cost function can be made ~0 and choose t based on this
    % distance.
    t=mu*M/(norm(A*x-b)+la*norm(x,1));
    ep=1e-10;
    al=0.15;
    be=0.8;

    fopts{1}=t;

    for it=1:maxiter
        z=bpdn_ncs(z,ep,al,be,t,maxiter);
        if (DEBUG)
            x=z(1:N);
            display(sprintf('outer iteration: %d',it));
            display(sprintf('Current cost: %f',norm(A*x-b)+la*norm(x,1)));
        end
        if M/t < ep
            break
        end
        t=mu*t;
        fopts{1}=t;
    end

    function [p]=f(z,fopts)
        % Evaluate cost function
        t=fopts{1};
        x=z(1:N);
        u=z(N+1:2*N);
        p=A*x;
        p=t*(x'*A'*p-2*b'*p+b'*b+la*sum(u));
        p=p-sum(log(u-x))-sum(log(x+u))-sum(log(u));
    end

    function [s]=bpdn_btls(z,dz,df,al,be,fopts,maxiter=100)
        % put z+s*dz in the domain of f
        s=1; 
        for it=1:maxiter
            z_=z+s*dz;
            x_=z_(1:N);
            u_=z_(N+1:2*N);
            if all((u_-x_)>0) && all((u_+x_)>0) && all(u_>0)
                break;
            end
            s*=be;
        end
        p=f(z,fopts);
        while (f(z+s*dz,fopts)>(p+al*s*df'*dz)) && (s > 0.01)
            s=be*s;
        end
        if s < 0.00001
            error('Step size is 0.');
        end
    end

    function [z]=bpdn_ncs(z,ep,al,be,t,maxiter)
        % Newton centering step
        fopts{1}=t;
        for it=1:maxiter
            x=z(1:N);
            u=z(N+1:2*N);
            d1_=u-x;
            d2_=u+x;
            d3_=u;
            d1=d1_.^-1;
            d2=d2_.^-1;
            d3=d3_.^-1;
            D1=diag(d1_.^-2);
            D2=diag(d2_.^-2);
            D3=diag(d3_.^-2);
            D=(1/(2*t))*(D1+D2-(D2-D1)*(D1+D2+D3)^(-1)*(D2-D1));
            rhs=-t*(2*A'*(A*x)-2*A'*b);
            rhs=rhs-d1+d2+(D2-D1)*(D1+D2+D3)^(-1)*(t*la*ones(N,1)-d1-d2-d3);
            r=rhs/(2*t);
            A_=D^(-1/2)*A';
            S=eye(M)+A_'*A_;
            R=chol(S);
            v=linsolve(R,...
                linsolve(R',A*D^(-1)*r,struct('LT','true')),...
                struct('UT','true'));
            dx=D^(-1)*(r-A'*v);
            du=(D1+D2+D3)^(-1)*(-t*la*ones(N,1)+d1+d2+d3-(D2-D1)*dx);
            gf=[t*(2*A'*(A*x)-2*A'*b)+d1-d2;t*la*ones(N,1)-d1-d2-d3];
            %d1_=u-x;
            %d2_=u+x;
            %d1=d1_.^-1;
            %d2=d2_.^-1;
            %D1=diag(d1_.^-2);
            %D2=diag(d2_.^-2);
            %g1=A*x;
            %g1=A'*g1-A'*b;
            %g2=t*la*ones(N,1)-d1-d2;
            %D=(D1+D2)^(-1)*2*D1*D2/t;
            %r=-g1+(1/(t*2))*(d2-d1+(D2-D1)*(D1+D2)^(-1)*g2);
            %v=linsolve(A*D^(-1)*A'+eye(M),A*D^(-1)*r);
            %dx=D^(-1)*(r-A'*v);
            %du=(D1+D2)^(-1)*(-g2-(D2-D1)*dx);
            dz=[dx;du];
            %gf=[2*t*g1+d1-d2;g2];
            lamb2=-gf'*dz;
            if (DEBUG)
                display(sprintf('inner iteration: %d',it));
                display(sprintf('lambda^2: %f',lamb2));
            end
            if (lamb2/2<=ep)
                break;
            end
            s=bpdn_btls(z,dz,gf,al,be,fopts);
            if (DEBUG)
                display(sprintf('line search result: %f',s));
            end
            z=z+s*dz;
        end
    end
end
