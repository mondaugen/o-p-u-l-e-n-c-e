function []=rand_par_cxn_file(N=100,...
                              pmin=5,...
                              pmax=10,...
                              path='/tmp/rand_%s.raw',...
                              dtime=1)
% RAND_PAR_CXN_FILE
%
% Generate random partials file and random connection file
%
% N is number of connection vectors
% pmin is minimum number of partials at each time step 
% pmax is maximum number of partials at each time step
% path is a string that contains one %s, which will be filled in with
%   {par,cxn} and will be used as the path, e.g., '/tmp/rand_%s.raw' will
%   generate the paths '/tmp/rand_par.raw' and '/tmp/rand_cxn.raw'
% dtime is the amount of time in seconds between partial vectors
%
L_PARTIAL_RECORD=4;
data=[];
lastData=[];
f=fopen(sprintf(path,'par'),'w');
g=fopen(sprintf(path,'cxn'),'w');
for k=1:N
    p=randi([pmin,pmax],1);
    % First two are amp, freq, rest are damp, dfreq
    data=[rand(p,1) sort(rand(p,1)) randn(p,2)];
    fwrite(f,dtime,'double');
    len_=p*L_PARTIAL_RECORD;
    fwrite(f,len_,'uint32');
    fwrite(f,data(:),'double');
    if (length(lastData)==0)
        lastData=data;
        continue;
    end
    m=size(data,1);
    n=size(lastData,1);
    x=rand(n*m,1);
    fwrite(g,m,'uint32');
    fwrite(g,n,'uint32');
    fwrite(g,x,'double');
    lastData=data;
end
fclose(f);
fclose(g);
