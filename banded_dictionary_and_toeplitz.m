ev=[(0:(N/8-1))/(N/8) 2-N/8:(N/4-1)/(N/8)]; 
plot(1:length(ev),ev);
ev=[(0:(N/8-1))/(N/8) 2-(N/8:(N/4-1))/(N/8)]; 
plot(1:length(ev),ev);
ev=[(0:(N/8-1))/(N/8) 2-(N/8:(N/4-1))/(N/8) zeros(1,N-N/4)]; 
octave:91> length(ev)
ans =  256
octave:92> plot(1:length(ev),ev);
octave:93> y=x.*ev;
octave:94> plot(1:length(y),y);
octave:95> y(64)
ans =  0.064145
octave:96> y(65)
ans = 0
octave:97> y_=circshift(y,[0 256-64]);
octave:98> plot(1:length(y_),y_);
octave:99> Y=zeros(256-64+1,length(y));
octave:100> K=(256-64);
octave:101> for k=0:K
> Y(k+1,:)=circshift(y,[0 k]);
