% Script file : bssprocess8.m

%-------------------------------------------------------
%  Frequency Domain Blind Source Separation 
%                  date: 12/10 2007, Masahiro Onohara
%
%  [2-input 2-output] freqency domain blind source separation.
%
%  Reference
%  Murata et al. (2001): "An Approach to Blind Source Separation Based on Temporal Structure of Speech Signalsz"
%
%  input parameters, there are 3.
%		x1	a observed signal (mixing signal)
%		x2	the other observed signal (mixing signal)
%		ilist	[fs; NFFT; shiftT; learntime]
%			fs:	sampling freqency 
%			NFFT:	length NFFT discrete Fourier transforms
%			shiftT:	window shift size
%			learntime:	learning time for unmixing [sec]
%
%  output parameters
%		y1a,y2a	and y1b,y2b	two pairs of unmixing signal, time domain	
%		Y1a,Y2a and Y1b,Y2b	two pairs of unmixing signal, freqency domain




% input section
% input list: fs, x1, x2, learntime, NFFT, shiftT, ilist

bssinput

%------------------

% STFT section

disp('STFT section');
tic

%---------
% STFT

X1=stft(x1,NFFT,fs,shiftT);
X2=stft(x2,NFFT,fs,shiftT);

sample_number=size(X1,2);	        % STFT後のサンプル数

if learntime==0, learn_number=0;			% learn_numberは学習サンプル数
else	learn_number=floor(2*learntime*fs/NFFT)-1;	% ０のときは全データ学習に用いるようにしている
end

% 直流成分bin=1は除去(=0)
% 以降、直流成分bin=1は除いて考える。
X1(1,:)=0;
X2(1,:)=0;

% NFFT/2以降の対象部分をカットして高速化。最後に復元。
X1(NFFT/2+2:NFFT,:)=[];
X2(NFFT/2+2:NFFT,:)=[];

% zero mean
for bin=2:NFFT/2+1
	X1(bin,:)=X1(bin,:)-mean(X1(bin,:));
	X2(bin,:)=X2(bin,:)-mean(X2(bin,:));	
end

disp('clear');
toc

%------------------

% learning section

disp('learning section');
tic

% iteration initialize
U1=zeros(size(X1));			% U1,U2は、分離されたスペクトルデータの一時的な保管場所
U2=zeros(size(X2));
B1=zeros(size(X1,1),2);	              	% Bは、求められた分離ベクトルwvの保管場所
B2=zeros(size(X2,1),2);

% iteration

for bin=2:NFFT/2+1

	% learning data
	if learn_number==0,				% learn_number=0のときは、全データで学習
		m1=X1(bin,:).';
		m2=X2(bin,:).';
	else	m1=X1(bin,1:learn_number).';		% 横ベクトルを縦ベクトルへ転置
		m2=X2(bin,1:learn_number).';		% mは学習用データ
	end

	%B=unmixing(m1,m2);	      			% 分離行列B=WV returned.
	B=timebss(m1,m2,40);

	u1=B(1,1)*X1(bin,:)+B(1,2)*X2(bin,:);		% uは分離信号	
	u2=B(2,1)*X1(bin,:)+B(2,2)*X2(bin,:);		% しかし置換と大きさの任意性あり

	U1(bin,:) = [u1];		% 分離信号を横ベクトルで収納
	U2(bin,:) = [u2];
	
	B1(bin,:) = B(1,:);		% 分離ベクトルを収納
	B2(bin,:) = B(2,:);

end

clear m1 m2 u1 u2 B bin learn_number

disp('clear');
toc

%----------------

% permutation and scaling section

disp('permutation section');
tic

[Y1a,Y1b,Y2a,Y2b] = permutation(U1,U2,B1,B2,ilist);

test2(Y1a,Y2a,NFFT)

clear U1 U2 B1 B2

disp('clear');
toc

%---------------------------------------
%---------------------------------------

% iSTFT section

% NFFT/2+1以降の周波数を復元
for m=0:NFFT/2-2,
	Y1a=[Y1a;Y1a(NFFT/2-m,:)];
	Y2a=[Y2a;Y2a(NFFT/2-m,:)];
	Y1b=[Y1b;Y1b(NFFT/2-m,:)];
	Y2b=[Y2b;Y2b(NFFT/2-m,:)];
end

disp('iSTFT section');
tic

y1a=istft(Y1a,NFFT,fs,shiftT);
y2a=istft(Y2a,NFFT,fs,shiftT);
y1b=istft(Y1b,NFFT,fs,shiftT);
y2b=istft(Y2b,NFFT,fs,shiftT);

disp('clear');
toc

x1=x1./max(abs(x1));
x2=x2./max(abs(x2));
y1a=y1a./max(abs(y1a));
y2a=y2a./max(abs(y2a));
y1b=y1b./max(abs(y1b));
y2b=y2b./max(abs(y2b));

figure;
subplot(2,3,1); plot(x1); xlabel('time'); ylabel('amp'); title('x1 signal');
subplot(2,3,2); plot(y1a); xlabel('time'); ylabel('amp'); title('y1a signal');
subplot(2,3,3); plot(y1b); xlabel('time'); ylabel('amp'); title('y1b signal');
subplot(2,3,4); plot(x2); xlabel('time'); ylabel('amp'); title('x2 signal');
subplot(2,3,5); plot(y2a); xlabel('time'); ylabel('amp'); title('y2a signal');
subplot(2,3,6); plot(y2b); xlabel('time'); ylabel('amp'); title('y2b signal');

figure;
subplot(2,3,1); specgram(x1,1024,fs); xlabel('time'); ylabel('freqency'); title('x1 signal');
subplot(2,3,2); specgram(y1a,1024,fs); xlabel('time'); ylabel('freqency'); title('y1a signal');
subplot(2,3,3); specgram(y1b,1024,fs); xlabel('time'); ylabel('freqency'); title('y1b signal');
subplot(2,3,4); specgram(x2,1024,fs); xlabel('time'); ylabel('freqency'); title('x2 signal');
subplot(2,3,5); specgram(y2a,1024,fs); xlabel('time'); ylabel('freqency'); title('y2a signal');
subplot(2,3,6); specgram(y2b,1024,fs); xlabel('time'); ylabel('freqency'); title('y2b signal');
	
clear fs NFFT learntime sample_number shiftT m
