function y = online_filtering(x,fs)

%-----------------------------------
%------------------------------------
% main program

% zero mean
x=x-mean(x);

% 基本周波数の計算（計算のシフト幅は10msec）
[t,f0]=shrp(x,fs);		% f0に区間毎の基本周波数の情報

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
for k=1:length(t),
	p=t(k);
	T(p:p+9)=1;
end

clear k p sa t va
% returned list
% x(対象信号)
% fs(サンプリング周波数)
% f0(区間毎の基本周波数情報、シフト幅10msec、計算範囲40msec)
% T(矩形信号 T(k)=1:ピッチ一定でフィルタリング要。T(k)=0:フィルタリング不要) 

%-----------------------
% filtering
int=0.010*fs/2;			% shift_interval シフト幅/2
first_shift=0.040*fs/2;		% 始めの計算範囲分をシフト（40msecの半分）
%y=zeros(size(x));		% フィルタリング後の信号を収納

% 最初の区間
t=1;
if T(t)==0, y(1:first_shift+int)=x(1:first_shift+int);	% フィルタリングなし
else 	ff=f0(t);					% フィルタリングあり、その区間の基本周波数
	y(1:first_shift+int)=zero_phase_filt(x(1:first_shift+int),ff,fs);	% フィルタリング
end

% 2番目以降の区間
for t=2:length(f0)-1,
	ts=2*(t+1)*int;		% 区間の中心の時点、この前後shiftの範囲を対象とする
	if T(t)==0, y(ts-int+1:ts+int)=x(ts-int+1:ts+int);		% フィルタリングなし
	else	ff=f0(t);			% その区間の基本周波数
		y(ts-int+1:ts+int)=zero_phase_filt(x(ts-int+1:ts+int),ff,fs);	% フィルタリング
	end
end

% 最後の区間
t=t+1;
ts=2*(t+1)*int;		% 区間の中心の時点、この前後shiftの範囲を対象とする
if T(t)==0, y(ts-int+1:length(x))=x(ts-int+1:end);		% フィルタリングなし
else	ff=f0(t);			% その区間の基本周波数
	y(ts-int+1:length(x))=zero_phase_filt(x(ts-int+1:end),ff,fs);	% フィルタリング
end

% normalize
y=y-mean(y);
y=y/max(abs(y));
y=y';

t=1:length(x);t=t./fs;
figure;
subplot(3,2,1);specgram(x,1024,16000,hann(1024),1000)
subplot(3,2,2);specgram(y,1024,16000,hann(1024),1000)
subplot(3,2,3);bar(T);axis([0,length(T),0,1]);
subplot(3,2,4);bar(T);axis([0,length(T),0,1]);
subplot(3,2,5);plot(t,x);axis([0,length(x)/fs,-1,1])
subplot(3,2,6);plot(t,y);axis([0,length(x)/fs,-1,1])

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
