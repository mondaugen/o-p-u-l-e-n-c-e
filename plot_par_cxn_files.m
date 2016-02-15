function [cxns]=plot_par_cxn_files(par_path,cxn_path,thresh=0.1,fig_num=1);
% PLOT_PAR_CXN_FILES
%
% par_path - path to partials file
% cxn_path - path to connections file
% thresh   - if a connection weight is over this value, plot it.
%
L_PARTIAL_RECORD=4;
f=fopen(par_path,'r');
g=fopen(cxn_path,'r');
newplot(figure(fig_num));
figure(fig_num);
hold on;
data=[];
lastData=[];
time=0;
cxns={};
cxns_idx=1;
while ~(feof(f) | feof(g))
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
    if (length(lastData)==0)
        lastData=data;
        continue;
    end
    m=fread(g,1,'uint32');
    n=fread(g,1,'uint32');
    cxn=fread(g,m*n,'double');
    cxn=reshape(cxn,[n,m])';
    cxns{cxns_idx}=cxn;
    cxns_idx+=1;
    plot(time,lastData(:,2));
    time+=dtime;
    plot(time,data(:,2));
    m;
    n;
    size(lastData,1);
    size(data,1);
    for i=(1:m)
        for j=(1:n)
            if (cxn(i,j) > thresh)
                % plot s.t. y value is freq 
                plot([time-dtime,time],[lastData(j,2),data(i,2)]);
            end
        end
    end
    lastData=data;
end
hold off;
