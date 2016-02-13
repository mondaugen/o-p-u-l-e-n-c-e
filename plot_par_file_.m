function []=plot_par_file(par_path,thresh=0.1);
% PLOT_PAR_CXN_FILES
%
% par_path - path to partials file
% cxn_path - path to connections file
% thresh   - if a connection weight is over this value, plot it.
%
L_PARTIAL_RECORD=4;
f=fopen(par_path,'r');
figure 1;
hold on;
data=[];
lastData=[];
time=0;
while ~(feof(f))
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
    plot(time,lastData(:,2));
    time+=dtime;
    plot(time,data(:,2));
    lastData=data;
end
hold off;
