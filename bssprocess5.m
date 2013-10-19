% Script file : bssprocess5.m

%---------------------------------
%---------------------------------

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

U=eigenvalue_decomposition(W*W');	% U = W*W'^(-1/2) = E*(D^(-1/2))*E';

W=U*W;	%orthogonalization

w1=W(1,:)';
w2=W(2,:)';

% learning

% iteration initialize

counter=0;

y1=w1(1)*m1+w1(2)*m2;
y2=w2(1)*m1+w2(2)*m2;
l=length(y1);
u=randn(l,1);			% gauss distribution u
Hu=-exp((-u.^2)/2);
neg=negentropy(y1,y2,l,Hu);

figure(1); hold on; plot(counter,neg,'.');
xlabel('counter'); ylabel('negentropy'); title('negentropy convergence');

%iteration

while counter<30 
	
	counter=counter+1;

	w1=bsslearning(m1,m2,w1);
	w2=bsslearning(m1,m2,w2);

	W=[w1 w2]';

	U=eigenvalue_decomposition(W*W');
	W=U*W;

	w1=W(1,:)';
	w2=W(2,:)';
	
	y1=w1(1)*m1+w1(2)*m2;
	y2=w2(1)*m1+w2(2)*m2;
	neg=negentropy(y1,y2,l,Hu);
	plot(counter,neg,'.');

end

hold off;

disp('unmixing matrix W*V');
W*V

clear counter m1 m2 l Hu u y1 y2 neg U

% unmixing signal y

y1=w1(1)*z1+w1(2)*z2;
y2=w2(1)*z1+w2(2)*z2;

%----------------------------------
%----------------------------------
%figure

x1=x1/(var(x1))^(1/2);
x2=x2/(var(x2))^(1/2);

figure(2)
subplot(1,3,1); plot(x1(1:8:16000),x2(1:8:16000),'.');
axis('square'); xlabel('x1'); ylabel('x2'); title('x1-x2');
subplot(1,3,2); plot(z1(1:8:16000),z2(1:8:16000),'.');
axis('square'); xlabel('z1'); ylabel('z2'); title('z1-z2');
subplot(1,3,3); plot(y1(1:8:16000),y2(1:8:16000),'.');
axis('square'); xlabel('y1'); ylabel('y2'); title('y1-y2');

x1=x1./max(abs(x1));
x2=x2./max(abs(x2));
y1=y1./max(abs(y1));
y2=y2./max(abs(y2));

figure(3)
subplot(2,2,1); plot(x1); xlabel('time'); ylabel('x1'); title('x1 signal');
subplot(2,2,2); plot(y1); xlabel('time'); ylabel('y1'); title('y1 signal');
subplot(2,2,3); plot(x2); xlabel('time'); ylabel('x2'); title('x2 signal');
subplot(2,2,4); plot(y2); xlabel('time'); ylabel('y2'); title('y2 signal');

%------------------------------

clear z1 z2 w1 w2 lim

disp('original_signal s1 s2, observed_signal x1 x2, unmixing_signal y1 y2, A,V,W reterned.');
