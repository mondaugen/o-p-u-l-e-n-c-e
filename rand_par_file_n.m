function []=rand_par_file_n(N=100,...
                              p=10,...
                              path='/tmp/rand_%s_n.raw',...
                              dtime=1,
                              sigma=0.01)
% RAND_PAR_CXN_FILE_N
%
% Generate random partials file whose partials do a random walk in frequency
%
% N is number of connection vectors so N+1 partial vectors will be produced
% p is number of partials at each time step 
% path is a string that contains one %s, which will be filled in with
%   {par,cxn} and will be used as the path, e.g., '/tmp/rand_%s.raw' will
%   generate the path '/tmp/rand_par.raw'
% dtime is the amount of time in seconds between partial vectors
% sigma is standard variation for the gaussian distributing the amount by which
% the frequencies change each time step
%
L_PARTIAL_RECORD=4;
data=[];
lastData=[];
f=fopen(sprintf(path,'par'),'w');
start_freqs=sort(rand(p,1));
for k=1:N
    % First two are amp, freq, rest are damp, dfreq
    dfreq=randn(p,1)*sigma;
    data=[rand(p,1) start_freqs randn(p,1) dfreq];
    start_freqs+=dfreq;
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
    lastData=data;
end
fclose(f);
