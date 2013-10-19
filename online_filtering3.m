function [y1,y2] = online_filtering(x1,x2,fs)

%-----------------------------------
%------------------------------------
% main program

% zero mean
x1=x1-mean(x1);
x2=x2-mean(x2);

% 基本周波数の計算（計算のシフト幅は10msec）
[t,f1]=shrp(x1,fs);		% f0に区間毎の基本周波数の情報
[t,f2]=shrp(x2,fs);
f0=(f1+f2)/2;

%--------------------------
% 一定ピッチの存在区間を求める

% 基本周波数の差（変化率）計算
sa=[0];
for k=2:length(f0),
	sa(k)=f0(k)-f0(k-1);
end

% 分散の計算範囲を100msecで計算。つまりf0で10点分。
for k=1:length(f0)-10,
	va(k)=var(sa(k:k+10));
end

% 矩形信号を作成（ピッチ一定:1,一定でない:0）
for k=1:length(f0)-10,
	if va(k)>200, T(k)=0;
	else T(k)=1;
	end
end
T=[T zeros(1,10)];
t=[];
for k=2:length(f0),
	if T(k)<T(k-1), t=[t;k];
	end
end
for k=1:length(t)-1,
	p=t(k);
	T(p:p+10)=1;
end

clear k p sa t va

% ピッチ一定区間の始まりと終わりの点を調べる
pre=0;
D=[];
for k=1:length(T),
	if sum([pre T(k)]==[0 1])==2, D=[D;k];
	elseif sum([pre T(k)]==[1 0])==2, D=[D;k];
	end
	pre=T(k);
end
if T(length(T))==1, D=[D;length(T)]; end
D=reshape(D,2,length(D)/2)';

% 一定区間の近似基本周波数を求める。区間内の平均と言う意味で近似。
F=[];
for k=1:size(D,1),
	F=[F; mean(f0(D(k,1):D(k,2)))];
end

% returned list
% x(対象信号)
% fs(サンプリング周波数)
% f0(区間毎の基本周波数情報、シフト幅10msec、計算範囲40msec)
% T(矩形信号 T(k)=1:ピッチ一定でフィルタリング要。T(k)=0:フィルタリング不要) 

%-----------------------
% filtering
int=0.010*fs/2;			% shift_interval シフト幅(10msec)/2
first_shift=0.040*fs/2;		% 始めの計算範囲分をシフト（40msecの半分）
y1=x1;
y2=x1;

for t=1:length(F),
	t1=int*(2*D(t,1)+1)+1;
	t2=int*(2*D(t,2)+3);
	y1(t1:t2)=zero_phase_filt(x1(t1:t2),F(t),fs);	% フィルタリング
	y2(t1:t2)=zero_phase_filt(x2(t1:t2),F(t),fs);
end

% normalize
y1=y1-mean(y1);
y2=y2-mean(y2);
y1=y1/max(abs(y1));
y2=y2/max(abs(y2));

%t=1:length(x1);t=t./fs;
%figure;
%subplot(3,2,1);specgram(x1,1024,16000,hann(1024),1000)
%subplot(3,2,2);specgram(y1,1024,16000,hann(1024),1000)
%subplot(3,2,3);bar(T);axis([0,length(T),0,1]);
%subplot(3,2,4);bar(T);axis([0,length(T),0,1]);
%subplot(3,2,5);plot(t,x1);axis([0,length(x1)/fs,-1,1])
%subplot(3,2,6);plot(t,y1);axis([0,length(x1)/fs,-1,1])
%figure;
%subplot(3,2,1);specgram(x2,1024,16000,hann(1024),1000)
%subplot(3,2,2);specgram(y2,1024,16000,hann(1024),1000)
%subplot(3,2,3);bar(T);axis([0,length(T),0,1]);
%subplot(3,2,4);bar(T);axis([0,length(T),0,1]);
%subplot(3,2,5);plot(t,x2);axis([0,length(x2)/fs,-1,1])
%subplot(3,2,6);plot(t,y2);axis([0,length(x2)/fs,-1,1])

%------------------------------------
%------------------------------------
% make filter and filtering
function y=zero_phase_filt(x,ff,fs)

%zero
zz=[];
ampz=[0.999 0.999 0.999 0.99];
for k=1:4,
	zz=[zz; ampz(k)*exp(i*2*pi*k*ff/fs); ampz(k)*exp(-i*2*pi*k*ff/fs)];
end

%pole
zp=[];
ampp=[0.92 0.985 0.985 0.99];
for k=1:4,
	zp=[zp; ampp(k)*exp(i*2*pi*(k-0.5)*ff/fs); ampp(k)*exp(-i*2*pi*(k-0.5)*ff/fs)];
end

b=poly(zz);	% zero
a=poly(zp);	% pole
y=filter(b,a,x);
if max(abs(y))>1, y=y./max(abs(y)); end
