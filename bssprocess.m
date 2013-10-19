% Script file : bssprocess.m

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

cov(z1,z2)

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

while counter<10 
	
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

disp('Evaluation Measure: W*V*A (normalized)');
WVA=W*V*A;
WVA(1,:)=WVA(1,:)/norm(WVA(1,:));
WVA(2,:)=WVA(2,:)/norm(WVA(2,:));
WVA

clear counter m1 m2 l Hu u y1 y2 neg U

% unmixing signal y

y1=w1(1)*z1+w1(2)*z2;
y2=w2(1)*z1+w2(2)*z2;

%----------------------------------
%----------------------------------
%figure

s1=s1/(var(s1))^(1/2);
s2=s2/(var(s2))^(1/2);

x1=x1/(var(x1))^(1/2);
x2=x2/(var(x2))^(1/2);

figure(2)
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

figure(3)
subplot(2,3,1); plot(s1); xlabel('time'); ylabel('s1'); title('s1 signal');
subplot(2,3,2); plot(x1); xlabel('time'); ylabel('x1'); title('x1 signal');
subplot(2,3,3); plot(y1); xlabel('time'); ylabel('y1'); title('y1 signal');
subplot(2,3,4); plot(s2); xlabel('time'); ylabel('s2'); title('s2 signal');
subplot(2,3,5); plot(x2); xlabel('time'); ylabel('x2'); title('x2 signal');
subplot(2,3,6); plot(y2); xlabel('time'); ylabel('y2'); title('y2 signal');

%------------------------------

clear z1 z2 w1 w2 lim

disp('original_signal s1 s2, observed_signal x1 x2, unmixing_signal y1 y2, A,V,W reterned.');
