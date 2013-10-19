% Script file : bssprocess2.m

%----------------------------------
%----------------------------------

fs=input('sampling freqency fs: ');
x1=input('Enter a original signal x1: ');
x2=input('Enter the other original signal x2: ');
learntime=input('learning and unmixing time [sec] : ');

x1=100*x1;
x2=100*x2;

% whitening

tic

% m is learning data for whitening
m1=x1(1:learntime*fs);
m2=x2(1:learntime*fs);

% zero mean
m1=m1-mean(m1);
m2=m2-mean(m2);

% decorrelation
Cov=cov(m1,m2);	
V=eigenvalue_decomposition(Cov);

% decorrelation
z1=V(1,1)*x1+V(1,2)*x2;
z2=V(2,1)*x1+V(2,2)*x2; 

% normalize var(z)=1
z1=z1/(var(z1))^(1/2);
z2=z2/(var(z2))^(1/2);

clear Cov m1 m2

%--------------------------------------
%--------------------------------------
% unmixing

% m is the learning date for unmixing.
m1=z1(1:learntime*fs); m2=z2(1:learntime*fs);

% zero mean
m1=m1-mean(m1);
m2=m2-mean(m2);

% initialize
w1=randn(2,1);
w2=randn(2,1);
w1=w1/norm(w1);
w2=w2/norm(w2);
W=[w1 w2]';

% Orthogonalization, eigenvalue decomposition W*W'
U=eigenvalue_decomposition(W*W');
W=U*W;

% learning

% iteration initialize

count=0;

y1=W(1,1)*m1+W(1,2)*m2;
y2=W(2,1)*m1+W(2,2)*m2;
l=length(y1);
WW=zeros(2,2);

%iteration

while log10(max(abs(W-WW)))>-3

	if count==30,break,end

	count=count+1;
	WW=W;

	gy1=tanh(y1);
	gy2=tanh(y2);

	beta1=-(y1'*gy1)/l;
	beta2=-(y2'*gy2)/l;

	arfa1=-1/( beta1 + sum(1./((cosh(y1)).^2))/l );
	arfa2=-1/( beta2 + sum(1./((cosh(y2)).^2))/l );
	
	R=([gy1'*y1 gy1'*y2; gy2'*y1 gy2'*y2])./l;
	
	W = W + (diag([arfa1 arfa2]))*([diag([beta1 beta2]) + R])*(W);

	U=eigenvalue_decomposition(W*W');
	W=U*W;

	y1=W(1,1)*m1+W(1,2)*m2;
	y2=W(2,1)*m1+W(2,2)*m2;
	
end

toc

count

disp('unmixing matrix W*V');
WV=W*V;
WV

clear arfa1 arfa2 beta1 beta2 gy1 gy2 m1 m2 y1 y2 U R l WW WV

% unmixing signal y

y1=W(1,1)*z1+W(1,2)*z2;
y2=W(2,1)*z1+W(2,2)*z2;

%----------------------------------
%----------------------------------
%figure

x1=x1/(var(x1))^(1/2);
x2=x2/(var(x2))^(1/2);

figure(3)
lim=16000;
subplot(1,3,1); plot(x1(1:8:lim),x2(1:8:lim),'.');
axis('square'); xlabel('x1'); ylabel('x2'); title('x1-x2');
subplot(1,3,2); plot(z1(1:8:lim),z2(1:8:lim),'.');
axis('square'); xlabel('z1'); ylabel('z2'); title('z1-z2');
subplot(1,3,3); plot(y1(1:8:lim),y2(1:8:lim),'.');
axis('square'); xlabel('y1'); ylabel('y2'); title('y1-y2');

x1=x1./max(abs(x1));
x2=x2./max(abs(x2));
y1=y1./max(abs(y1));
y2=y2./max(abs(y2));

figure(4)
subplot(2,2,1); plot(x1); xlabel('time'); ylabel('x1'); title('x1 signal');
subplot(2,2,2); plot(y1); xlabel('time'); ylabel('y1'); title('y1 signal');
subplot(2,2,3); plot(x2); xlabel('time'); ylabel('x2'); title('x2 signal');
subplot(2,2,4); plot(y2); xlabel('time'); ylabel('y2'); title('y2 signal');

%------------------------------

clear z1 z2 w1 w2 lim W11 W12 W21 W22 WVA11 WVA12 WVA21 WVA22 itercount count 

disp('original_signal s1 s2, observed_signal x1 x2, unmixing_signal y1 y2, A,V,W reterned.');

