% Script file : bssprocess8.m

%-------------------------------------
%-------------------------------------
% input section
% input list: fs, x1, x2, learntime, NFFT, shiftT, ilist
bssinput

%------------------
% STFT section

disp('STFT section');
tic

X1=stft(x1,NFFT,fs,shiftT);
X2=stft(x2,NFFT,fs,shiftT);

if learntime==0, T=size(X1,2);			% T�ϳؽ�����ץ��
else	T=floor((learntime*fs-NFFT)/shiftT+1);		% ���ΤȤ������ǡ����ؽ����Ѥ���褦�ˤ��Ƥ���
end

% NFFT/2�ʹߤ��о���ʬ�򥫥åȤ��ƹ�®�����Ǹ��������
X1(NFFT/2+2:NFFT,:)=[];
X2(NFFT/2+2:NFFT,:)=[];

disp('clear');
toc

%------------------
% learning section

disp('learning section');
tic

% iteration

for bin=1:NFFT/2+1,

	% learning data
	m1=X1(bin,1:T).';		% ���٥��ȥ��ĥ٥��ȥ��ž��
	m2=X2(bin,1:T).';		% m�ϳؽ��ѥǡ���

	%B=unmixing(m1,m2);	      			% ʬΥ����B=WV returned.
	B=timebss(m1,m2,50);
	u1=B(1,1)*X1(bin,:)+B(1,2)*X2(bin,:);		% u��ʬΥ����	
	u2=B(2,1)*X1(bin,:)+B(2,2)*X2(bin,:);		% �������ִ����礭����Ǥ��������

	% scaling�β��
	A=inv(B);					% ʬΥ����εչ��󡢤Ĥޤ꺮�������б�
	U1a(bin,:) = A(1,1)*u1;				% U1a��U2a��scaling���褷��ʬΥ����ΰ���
	U2a(bin,:) = A(1,2)*u2;				% U1b��U2b�Ǥ⤦����
	U1b(bin,:) = A(2,1)*u1;
	U2b(bin,:) = A(2,2)*u2;

end

clear m1 m2 u1 u2 B bin

disp('clear');
toc

%----------------

% permutation and scaling section

disp('permutation section');
tic

per = permutation2(U1a(:,1:T),U2a(:,1:T),U1b(:,1:T),U2b(:,1:T),ilist);

% initialize; permutation���褷��������Ǽ
Y1a=[];
Y2a=[];
Y1b=[];
Y2b=[];

% per�˽����ºݤ��¤��ؤ�
for bin=1:NFFT/2+1,
	if per(bin)==0,
		Y1a=[Y1a;U1a(bin,:)];
		Y2a=[Y2a;U2a(bin,:)];
		Y1b=[Y1b;U1b(bin,:)];
		Y2b=[Y2b;U2b(bin,:)];
	elseif	per(bin)==1,
		Y1a=[Y1a;U2a(bin,:)];
		Y2a=[Y2a;U1a(bin,:)];
		Y1b=[Y1b;U2b(bin,:)];
		Y2b=[Y2b;U1b(bin,:)];
	end
end

clear U1a U2a U1b U2b T per bin

disp('clear');
toc

%---------------------------------------
%---------------------------------------

% iSTFT section

disp('iSTFT section');
tic

% NFFT/2+1�ʹߤμ��ȿ�������
for m=0:NFFT/2-2,
	Y1a=[Y1a;Y1a(NFFT/2-m,:)];
	Y2a=[Y2a;Y2a(NFFT/2-m,:)];
	Y1b=[Y1b;Y1b(NFFT/2-m,:)];
	Y2b=[Y2b;Y2b(NFFT/2-m,:)];
end

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

tx=(1:length(x1))/fs;
ty=(1:length(y1a))/fs;

figure;
subplot(2,3,1); plot(tx,x1); xlabel('time[sec]'); ylabel('amp'); title('x1 signal');
subplot(2,3,2); plot(ty,y1a); xlabel('time[sec]'); ylabel('amp'); title('y1a signal');
subplot(2,3,3); plot(ty,y1b); xlabel('time[sec]'); ylabel('amp'); title('y1b signal');
subplot(2,3,4); plot(tx,x2); xlabel('time[sec]'); ylabel('amp'); title('x2 signal');
subplot(2,3,5); plot(ty,y2a); xlabel('time[sec]'); ylabel('amp'); title('y2a signal');
subplot(2,3,6); plot(ty,y2b); xlabel('time[sec]'); ylabel('amp'); title('y2b signal');

figure;
subplot(2,3,1); specgram(x1,1024,fs); xlabel('time'); ylabel('freqency'); title('x1 signal');
subplot(2,3,2); specgram(y1a,1024,fs); xlabel('time'); ylabel('freqency'); title('y1a signal');
subplot(2,3,3); specgram(y1b,1024,fs); xlabel('time'); ylabel('freqency'); title('y1b signal');
subplot(2,3,4); specgram(x2,1024,fs); xlabel('time'); ylabel('freqency'); title('x2 signal');
subplot(2,3,5); specgram(y2a,1024,fs); xlabel('time'); ylabel('freqency'); title('y2a signal');
subplot(2,3,6); specgram(y2b,1024,fs); xlabel('time'); ylabel('freqency'); title('y2b signal');

clear fs NFFT shiftT learntime tx ty m bin
