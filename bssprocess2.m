% Script file : bssprocess2.m

%-------------------------------
%-------------------------------
% mixing

fs=input('sampling freqency fs: ');
s1=input('Enter a original signal s1: ');
s2=input('Enter the other original signal s2: ');

% s times 10 for preventing calculation error
s1=10*s1;
s2=10*s2;

% make a same length signal
sa=length(s1)-length(s2);
if sa>0, s2 = [s2; zeros(sa,1)];
elseif sa<0, s1=[s1; zeros(-sa,1)];                     
end
clear sa;

% make mixing matrix
A=randn(2);

% make observed singnal x
x1=A(1,1)*s1+A(1,2)*s2;
x2=A(2,1)*s1+A(2,2)*s2;

disp(' ');
disp('mixing matrix A');
A

%----------------------------------
%----------------------------------
% whitening

tic

% m is learning data for whitening, use x's first 3 seconds
if length(x1)<3*fs, m1=x1; m2=x2;
else m1=x1(1:3*fs); m2=x2(1:3*fs);
end

% zero mean

m1=m1-mean(m1);
m2=m2-mean(m2);

% make decorrelation matrix

% eigenvalue decomposition 'Cov'

Cov=cov(m1,m2);	

% Cov=E*D*E'			eigenvalue decomposition
% V=E*(D^(-1/2))*E' 	decorrelation matrix

V=eigenvalue_decomposition(Cov);

% decorrelation
% get whitening data z 
% not learning data for whitening w but observed data x make z. assuming that x is zero mean.

z1=V(1,1)*x1+V(1,2)*x2;
z2=V(2,1)*x1+V(2,2)*x2; 

% normalize var(z)=1

z1=z1/(var(z1))^(1/2);
z2=z2/(var(z2))^(1/2);

clear Cov m1 m2

%--------------------------------------
%--------------------------------------
% unmixing

% unmixing vector w1 and w2.

% m is the learning date for unmixing.
if length(z1)<3*fs, m1=z1; m2=z2;
else m1=z1(1:3*fs); m2=z2(1:3*fs);
end

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

% learning

% iteration initialize

W11=[w1(1)];
W12=[w1(2)];
W21=[w2(1)];
W22=[w2(2)];

WVA=W*V*A;
WVA(1,:)=WVA(1,:)/norm(WVA(1,:));
WVA(2,:)=WVA(2,:)/norm(WVA(2,:));

WVA11=[WVA(1,1)];
WVA12=[WVA(1,2)];
WVA21=[WVA(2,1)];
WVA22=[WVA(2,2)];

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

	W11=[W11 W(1,1)];
	W12=[W12 W(1,2)];
	W21=[W21 W(2,1)];
	W22=[W22 W(2,2)];

	WVA=W*V*A;
	WVA(1,:)=WVA(1,:)/norm(WVA(1,:));
	WVA(2,:)=WVA(2,:)/norm(WVA(2,:));

	WVA11=[WVA11 WVA(1,1)];
	WVA12=[WVA12 WVA(1,2)];
	WVA21=[WVA21 WVA(2,1)];
	WVA22=[WVA22 WVA(2,2)];
	
	y1=W(1,1)*m1+W(1,2)*m2;
	y2=W(2,1)*m1+W(2,2)*m2;
	
end

toc

count

disp('unmixing matrix W*V');
WV=W*V;
WV

disp('Evaluation Measure: W*V*A (normalized)');
WVA=W*V*A;
WVA(1,:)=WVA(1,:)/norm(WVA(1,:));
WVA(2,:)=WVA(2,:)/norm(WVA(2,:));
WVA

clear arfa1 arfa2 beta1 beta2 gy1 gy2 m1 m2 y1 y2 U R l WW WV

% unmixing signal y

y1=W(1,1)*z1+W(1,2)*z2;
y2=W(2,1)*z1+W(2,2)*z2;

%----------------------------------
%----------------------------------
%figure

itercount=(0:1:count);

figure(1);
subplot(2,2,1); plot(itercount,W11,'.'); xlabel('itercount'); title('W(1,1)');
subplot(2,2,2); plot(itercount,W12,'.'); xlabel('itercount'); title('W(1,2)');
subplot(2,2,3); plot(itercount,W21,'.'); xlabel('itercount'); title('W(2,1)');
subplot(2,2,4); plot(itercount,W22,'.'); xlabel('itercount'); title('W(2,2)');
title('W convergence');

figure(2);
subplot(2,2,1); plot(itercount,WVA11,'.'); xlabel('itercount'); title('WVA(1,1)');
subplot(2,2,2); plot(itercount,WVA12,'.'); xlabel('itercount'); title('WVA(1,2)');
subplot(2,2,3); plot(itercount,WVA21,'.'); xlabel('itercount'); title('WVA(2,1)');
subplot(2,2,4); plot(itercount,WVA22,'.'); xlabel('itercount'); title('WVA(2,2)');
title('WVA convergence');

s1=s1/(var(s1))^(1/2);
s2=s2/(var(s2))^(1/2);

x1=x1/(var(x1))^(1/2);
x2=x2/(var(x2))^(1/2);

figure(3)
if length(x1)<16000, lim=length(x1);
else lim=16000;
end
subplot(2,2,1); plot(s1(1:8:lim),s2(1:8:lim),'.');
axis('square'); xlabel('s1'); ylabel('s2'); title('s1-s2');
subplot(2,2,2); plot(x1(1:8:lim),x2(1:8:lim),'.');
axis('square'); xlabel('x1'); ylabel('x2'); title('x1-x2');
subplot(2,2,3); plot(z1(1:8:lim),z2(1:8:lim),'.');
axis('square'); xlabel('z1'); ylabel('z2'); title('z1-z2');
subplot(2,2,4); plot(y1(1:8:lim),y2(1:8:lim),'.');
axis('square'); xlabel('y1'); ylabel('y2'); title('y1-y2');

s1=s1./max(abs(s1));
s2=s2./max(abs(s2));
x1=x1./max(abs(x1));
x2=x2./max(abs(x2));
y1=y1./max(abs(y1));
y2=y2./max(abs(y2));

figure(4)
subplot(2,3,1); plot(s1); xlabel('time'); ylabel('s1'); title('s1 signal');
subplot(2,3,2); plot(x1); xlabel('time'); ylabel('x1'); title('x1 signal');
subplot(2,3,3); plot(y1); xlabel('time'); ylabel('y1'); title('y1 signal');
subplot(2,3,4); plot(s2); xlabel('time'); ylabel('s2'); title('s2 signal');
subplot(2,3,5); plot(x2); xlabel('time'); ylabel('x2'); title('x2 signal');
subplot(2,3,6); plot(y2); xlabel('time'); ylabel('y2'); title('y2 signal');

%------------------------------

clear z1 z2 w1 w2 lim W11 W12 W21 W22 WVA11 WVA12 WVA21 WVA22 itercount count

disp('original_signal s1 s2, observed_signal x1 x2, unmixing_signal y1 y2, A,V,W reterned.');
