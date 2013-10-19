function [f0,sa,va]=timefilt(x,fs)

th=0.1;         % time threshold for constant pitch [sec]
timestep=0.01;  % timestep for analysis in pitch estimation
ts=fs*timestep; % 1sec = fs*timestep point

% pitch estimation
[a,f0]=shrp(x,fs);
%subplot(3,1,1);plot(f0);
clear a

sa=zeros(size(x));
for k=2:length(x),
    sa(k)=x(k)-x(k)-1;
end
%subplot(3,1,2);plot(sa);

for k=1:length(x)-ts*th,
    va(k)=var(x(k:k+ts*th));
end
%subplot(3,1,3);plot(va);


