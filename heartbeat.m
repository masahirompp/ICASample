function y=heartbeat(x)

% 前処理
x=x-mean(x);
x=x./max(abs(x));

y=abs(x);
m=sqrt(mean(y.^2))

% m以下はカット
for k=1:length(y),
	if y(k)<0.5, y(k)=0;	end
end
plot(y)
