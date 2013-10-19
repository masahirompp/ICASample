function[time,count,WVA,E1] = bssprocess3(fs,s1,s2,learntime)

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

%----------------------------------
%----------------------------------
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
m1=z1(1:learntime*fs);
m2=z2(1:learntime*fs);

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
w1=W(1,:)';
w2=W(2,:)';

% learning

% iteration initialize

count=0;
WW=zeros(2,2);

%iteration

while log10(max(abs(W-WW)))>-3

	if count==30,break,end
	
	count=count+1;
	WW=W;

	w1=bsslearning(m1,m2,w1);
	w2=bsslearning(m1,m2,w2);

	W=-[w1 w2]';

	U=eigenvalue_decomposition(W*W');
	W=U*W;

	w1=W(1,:)';
	w2=W(2,:)';

end

time=toc;
count;
WVA=W*V*A;
WVA(1,:)=WVA(1,:)/norm(WVA(1,:));
WVA(2,:)=WVA(2,:)/norm(WVA(2,:));
E1=evaluation(WVA);

clear m1 m2 l Hu u y1 y2 neg U w1 w2 z1 z2 WW