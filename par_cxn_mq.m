function [costs]=par_cxn_mq(par_path,cxn_path,mu)
    % PAR_CXN_MQ
    %
    % Solve partial connections using the technique of MacAulay and Quatieri
    % par_path - path to file containing partials
    % cxn_path - path to file where connections will be saved
    % mu       - threshold of max cost for path
    
    % Weight vector for amp,freq,damp,dfreq costs resp.
    L_PARTIAL_RECORD=4;
    w=[0;1;0;0];
    lastData=[];
    f=fopen(par_path,'r');
    g=fopen(cxn_path,'w');
    costs={};
    k=1;
    lastCosts=[];
    while ~feof(f)
        dtime=fread(f,1,'double');
        if(length(dtime)==0)
            break;
        end
        ndata=fread(f,1,'uint32');
        if(length(ndata)==0)
            break;
        end
        data=fread(f,ndata,'double');
        data=reshape(data,[length(data)/L_PARTIAL_RECORD,L_PARTIAL_RECORD]);
        % use only frequency data
        data=data(:,2);
        if (length(lastData)==0)
            lastData=data;
            continue;
        end
        m=size(data,1);
        n=size(lastData,1);
        C=pt_build_cost_matrix(lastData(:),data(:),@divergence);
        k+=1;
        % make sure mu is compatible with divergence
        xopt=pt_mq_cxns(C,divergence(mu,0));
        xopt=xopt';
        xopt=xopt(:);
        fwrite(g,m,'uint32');
        fwrite(g,n,'uint32');
        fwrite(g,xopt,'double');
        lastData=data;
    end
    fclose(f);
    fclose(g);
end
%
function [d]=divergence(x,y)
    d=(x.-y).^2;
end
