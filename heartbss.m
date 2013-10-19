function[x1,x2,y1,y2]=heartbss(ch1,ch2,fs,NFFT,shiftT,learntime);

%-------------------------------------
%-------------------------------------
% input section
% input list: fs, x1, x2, learntime, NFFT, shiftT

disp('filtering section');
tic

[x1,x2]=online_filtering3(ch1,ch2,fs);

X1=stft(x1,NFFT,fs,shiftT);
X2=stft(x2,NFFT,fs,shiftT);

if learntime==0, T=size(X1,2);			% Tは学習サンプル数
else	T=floor((learntime*fs-NFFT)/shiftT+1);		% ０のときは全データ学習に用いるようにしている
end

% NFFT/2以降の対象部分をカットして高速化。最後に復元。
X1(NFFT/2+2:NFFT,:)=[];
X2(NFFT/2+2:NFFT,:)=[];

disp('clear');
toc

%------------------
% learning section

disp('learning section');
tic

% スケーリングをch1とch2のどちらに合わせるか。ch1ならば1、ch2ならば2を入力。
scaling=2;

% iteration

for bin=1:NFFT/2+1,

	% learning data
	m1=X1(bin,1:T).';		% 横ベクトルを縦ベクトルへ転置
	m2=X2(bin,1:T).';		% mは学習用データ

	%B=unmixing(m1,m2);	      			% 分離行列B=WV returned.
	B=timebss(m1,m2,50);
	u1=B(1,1)*X1(bin,:)+B(1,2)*X2(bin,:);		% uは分離信号	
	u2=B(2,1)*X1(bin,:)+B(2,2)*X2(bin,:);		% しかし置換と大きさの任意性あり

	% scalingの解決
	A=inv(B);					% 分離行列の逆行列、つまり混合行列に対応
	U1(bin,:) = A(scaling,1)*u1;			% ch(scaling)にスケーリングを合わせる
	U2(bin,:) = A(scaling,2)*u2;

end

clear m1 m2 u1 u2 B bin

disp('clear');
toc

%----------------

% permutation and scaling section

disp('permutation section');
tic

per = permutation3(U1(:,1:T),U2(:,1:T),20);

% initialize; permutationを解決した信号を収納
Y1=[];
Y2=[];

% perに従い実際に並べ替え
for bin=1:NFFT/2+1,
	if per(bin)==0,
		Y1=[Y1;U1(bin,:)];
		Y2=[Y2;U2(bin,:)];
	elseif	per(bin)==1,
		Y1=[Y1;U2(bin,:)];
		Y2=[Y2;U1(bin,:)];
	end
end

clear U1 U2 T per bin

disp('clear');
toc

%---------------------------------------
%---------------------------------------

% iSTFT section

disp('iSTFT section');
tic

% NFFT/2+1以降の周波数を復元
for m=0:NFFT/2-2,
	Y1=[Y1;Y1(NFFT/2-m,:)];
	Y2=[Y2;Y2(NFFT/2-m,:)];
end

y1=istft(Y1,NFFT,fs,shiftT);
y2=istft(Y2,NFFT,fs,shiftT);

y1=y1./max(abs(y1));
y2=y2./max(abs(y2));

disp('clear');
toc

clear m bin Y1 Y2 M scaling

%---------------------
% figure section

tx=(1:length(x1))/fs;
ty=(1:length(y1))/fs;
figure;
subplot(2,3,1);plot(tx,ch1);axis([0,length(y1)/fs,-1,1]);title('ch1')
subplot(2,3,2);plot(tx,x1);axis([0,length(y1)/fs,-1,1]);title('x1')
subplot(2,3,3);plot(ty,y1);axis([0,length(y1)/fs,-1,1]);title('y1')
subplot(2,3,4);plot(tx,ch2);axis([0,length(y1)/fs,-1,1]);title('ch2')
subplot(2,3,5);plot(tx,x2);axis([0,length(y1)/fs,-1,1]);title('x2')
subplot(2,3,6);plot(ty,y2);axis([0,length(y1)/fs,-1,1]);title('y2')
figure;
subplot(2,3,1);specgram(ch1,1024,16000,hann(1024),1000);title('ch1')
subplot(2,3,2);specgram(x1,1024,16000,hann(1024),1000);title('x1')
subplot(2,3,3);specgram(y1,1024,16000,hann(1024),1000);title('y1')
subplot(2,3,4);specgram(ch2,1024,16000,hann(1024),1000);title('ch2')
subplot(2,3,5);specgram(x2,1024,16000,hann(1024),1000);title('x2')
subplot(2,3,6);specgram(y2,1024,16000,hann(1024),1000);title('y2')
