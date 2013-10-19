% Script file : bssprocess7.m

%-------------------------------------
%-------------------------------------

% input section
% input list: fs, x1, x2, learntime, NFFT, shiftT, ilist

bssinput

%------------------

% STFT section

disp('STFT section');
tic

%---------
% STFT

X1=specgram(x1,NFFT,fs,hann(NFFT),NFFT-shiftT);
X2=specgram(x2,NFFT,fs,hann(NFFT),NFFT-shiftT);

ans=size(X1);

bin_number=ans(1);		% ���ȿ�bin��
sample_number=ans(2);	        % STFT��Υ���ץ��

if learntime==0, learn_number=0;			% learn_number�ϳؽ�����ץ��
else	learn_number=floor(2*learntime*fs/NFFT)-1;	% ���ΤȤ������ǡ����ؽ����Ѥ���褦�ˤ��Ƥ���
end

% ľή��ʬbin=1�Ͻ���(=0)
% �ʹߡ�ľή��ʬbin=1�Ͻ����ƹͤ��롣
X1(1,:)=0;
X2(1,:)=0;

% zero mean
bin=1;
while bin < bin_number
	bin=bin+1;
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
bin=1;                    		% �׻���μ��ȿ�bin�ʥ�С���ľή��ʬ�Ͻ���
U1=zeros(size(X1));			% U1,U2�ϡ�ʬΥ���줿���ڥ��ȥ�ǡ����ΰ��Ū���ݴɾ��
U2=zeros(size(X2));
B1=zeros(bin_number,2);	              	% B�ϡ�����줿ʬΥ�٥��ȥ�wv���ݴɾ��
B2=zeros(bin_number,2);

testz=[0];
testy=[0];

% iteration

while bin < bin_number

	bin = bin+1;

	if learn_number==0,				% learn_number=0�ΤȤ��ϡ����ǡ����ǳؽ�
		m1=X1(bin,:).';
		m2=X2(bin,:).';
	else	m1=X1(bin,1:learn_number).';		% ���٥��ȥ��ĥ٥��ȥ��ž��
		m2=X2(bin,1:learn_number).';		% m�ϳؽ��ѥǡ���
	end

	[B,covz,covy]=unmixing(m1,m2);	      			% ʬΥ����B=WV returned.

	u1=B(1,1)*X1(bin,:)+B(1,2)*X2(bin,:);		% u��ʬΥ����	
	u2=B(2,1)*X1(bin,:)+B(2,2)*X2(bin,:);		% �������ִ����礭����Ǥ��������

	U1(bin,:) = [u1];		% ʬΥ����򲣥٥��ȥ�Ǽ�Ǽ
	U2(bin,:) = [u2];
	
	B1(bin,:) = B(1,:);		% ʬΥ�٥��ȥ���Ǽ
	B2(bin,:) = B(2,:);

	testz=[testz;covz];
	testy=[testy;covy];

end

figure;
subplot(2,1,1);plot(abs(testz))
subplot(2,1,2);plot(abs(testy))

clear m1 m2 u1 u2 B bin learn_number ans %X1 X2

disp('clear');
toc

%----------------

% permutation and scaling section

disp('permutation section');
tic

[Y1a,Y1b,Y2a,Y2b] = permutation(U1,U2,B1,B2,ilist);
%Y1a=U1;
%Y1b=U1;
%Y2a=U2;
%Y2b=U2;

%test2(Y1a,Y2a)

clear U1 U2 B1 B2

disp('clear');
toc

%---------------------------------------
%---------------------------------------

% iSTFT section

disp('iSTFT section');
tic

l=NFFT+shiftT*(sample_number-1);	% ���������Ĺ��
h=hann(NFFT);				% hanning window

%---------------
% hanning window �εմؿ���ȯ�����ʤ��褦����������������

H=zeros(l,1);
k=0;

while k<sample_number
	H(shiftT*k+1:shiftT*k+NFFT,1)=H(shiftT*k+1:shiftT*k+NFFT,1)+h;
	k=k+1;
end

H=(tanh(H-max(H))+max(H))./max(H);	% �����ʬΥ������뤳�Ȥǡ��뤫���εդ����򤹤뤳�Ȥˤʤ�

%-------------
% unmixing signal y1a, y2a, y1b, y2b 

y1a=zeros(l,1);
y2a=zeros(l,1);
y1b=zeros(l,1);
y2b=zeros(l,1);

ts=0;

while ts<sample_number

	m1a=real(ifft(Y1a(:,ts+1),NFFT));
	y1a(shiftT*ts+1:shiftT*ts+NFFT,1)=y1a(shiftT*ts+1:shiftT*ts+NFFT,1)+m1a;

	m2a=real(ifft(Y2a(:,ts+1),NFFT));
	y2a(shiftT*ts+1:shiftT*ts+NFFT,1)=y2a(shiftT*ts+1:shiftT*ts+NFFT,1)+m2a;

	m1b=real(ifft(Y1b(:,ts+1),NFFT));
	y1b(shiftT*ts+1:shiftT*ts+NFFT,1)=y1b(shiftT*ts+1:shiftT*ts+NFFT,1)+m1b;

	m2b=real(ifft(Y2b(:,ts+1),NFFT));
	y2b(shiftT*ts+1:shiftT*ts+NFFT,1)=y2b(shiftT*ts+1:shiftT*ts+NFFT,1)+m2b;

	ts=ts+1;

end

y1a=y1a./H;
y2a=y2a./H;
y1b=y1b./H;
y2b=y2b./H;

clear m1a m2a m1b m2b ts H h k l
%clear Y1a Y2a Y1b Y2b

disp('clear');
toc
	
clear fs NFFT learntime bin_number sample_number shiftT

x1=x1./max(abs(x1));
x2=x2./max(abs(x2));
y1a=y1a./max(abs(y1a));
y2a=y2a./max(abs(y2a));
y1b=y1b./max(abs(y1b));
y2b=y2b./max(abs(y2b));

figure;
subplot(2,3,1); plot(x1); xlabel('time'); ylabel('x1'); title('x1 signal');
subplot(2,3,2); plot(y1a); xlabel('time'); ylabel('y1a'); title('y1a signal');
subplot(2,3,3); plot(y1b); xlabel('time'); ylabel('y1b'); title('y1b signal');
subplot(2,3,4); plot(x2); xlabel('time'); ylabel('x2'); title('x2 signal');
subplot(2,3,5); plot(y2a); xlabel('time'); ylabel('y2a'); title('y2a signal');
subplot(2,3,6); plot(y2b); xlabel('time'); ylabel('y2b'); title('y2b signal');

