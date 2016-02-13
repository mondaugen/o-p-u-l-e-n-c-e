function []=rand_par_file_sin_noise(N=100,...
                              p=10,...
                              path='/tmp/rand_%s_sin.raw',...
                              dtime=1,
                              sigma=0.01,
                              fmodamp=0.1,
                              fmodfreq=0.1,
                              min_noise_p=0,
                              max_noise_p=10)
% RAND_PAR_CXN_FILE_SIN_NOISE
%
% Generate random partials file whose partials are modulated in frequency by a
% sine function and with some extra noise partials
%
% N is number of connection vectors so N+1 partial vectors will be produced
% p is number of partials at each time step 
% path is a string that contains one %s, which will be filled in with
%   {par,cxn} and will be used as the path, e.g., '/tmp/rand_%s.raw' will
%   generate the path '/tmp/rand_par.raw'
% dtime is the amount of time in seconds between partial vectors
% sigma is standard variation for the gaussian distributing the amount by which
% fmodamp is the amplitude of frequency modulation
% fmodfreq is the frequency of frequency modulation in periods/dtime
% the frequencies change each time step
% min_noise_p is the minimum number of noise partials produced at each time step
% max_noise_p is the maximum number of noise_partials produced at each time step
%
L_PARTIAL_RECORD=4;
data=[];
lastData=[];
f=fopen(sprintf(path,'par'),'w');
start_freqs=sort(rand(p,1));
time=0;
for k=1:N
    % First two are amp, freq, rest are damp, dfreq
    dfreq=ones(p,1)*fmodamp*fmodfreq*2*pi*dtime*cos(fmodfreq*2*pi*time);
    dfreq+=randn(p,1)*sigma;
    data=[rand(p,1) start_freqs randn(p,1) dfreq];
    noise_p=randi([min_noise_p max_noise_p],1);
    data=[data;
          rand(noise_p,2) randn(noise_p,2)];
    start_freqs+=dfreq;
    time+=dtime;
    fwrite(f,dtime,'double');
    len_=(p+noise_p)*L_PARTIAL_RECORD;
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
