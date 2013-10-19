%
%
% time dependent BSS algorythm

% Script file : timestrbss.m

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

% m is learning data for whitening, use x's first 2 seconds
if length(x1)<2*fs, m1=x1; m2=x2;
else m1=x1(1:2*fs); m2=x2(1:2*fs);
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

% m is the learning date for unmixing.
if length(z1)<3*fs, m1=z1; m2=z2;
else m1=z1(1:3*fs); m2=z2(1:3*fs);
end

% zero mean

m1=m1-mean(m1);
m2=m2-mean(m2);

l=length(m1);
tao=ceil(100*rand(1,1))		% ランダムな２桁の正の整数

%c11=([m1;zeros(tao,1)]'*[zeros(tao,1);m1])/(l-tao);
%c12=([m1;zeros(tao,1)]'*[zeros(tao,1);m2])/(l-tao);
%c21=([m2;zeros(tao,1)]'*[zeros(tao,1);m1])/(l-tao);
%c22=([m2;zeros(tao,1)]'*[zeros(tao,1);m2])/(l-tao);

c11=(m1(tao+1:l)'*m1(1:l-tao))/(l-tao);
c12=(m1(tao+1:l)'*m2(1:l-tao))/(l-tao);
c21=(m2(tao+1:l)'*m1(1:l-tao))/(l-tao);
c22=(m2(tao+1:l)'*m2(1:l-tao))/(l-tao);

C=[c11 c12; c21 c22]
C=(C+C')/2

[E,D]=eig(C);		% eigenvalue
W=E';

y1=W(1,1)*z1+W(1,2)*z2;
y2=W(2,1)*z1+W(2,2)*z2;


disp('unmixing matrix W*V');
W*V

disp('Evaluation Measure: W*V*A (normalized)');
WVA=W*V*A;
WVA(1,:)=WVA(1,:)/norm(WVA(1,:));
WVA(2,:)=WVA(2,:)/norm(WVA(2,:));
WVA


clear c11 c12 c21 c22 m1 m2 l tao C D E

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
axis([-8 8 -8 8]); xlabel('s1'); ylabel('s2'); title('s1-s2');
subplot(2,2,2); plot(x1(1:8:lim),x2(1:8:lim),'.');
axis([-8 8 -8 8]); xlabel('x1'); ylabel('x2'); title('x1-x2');
subplot(2,2,3); plot(z1(1:8:lim),z2(1:8:lim),'.');
axis([-8 8 -8 8]); xlabel('z1'); ylabel('z2'); title('z1-z2');
subplot(2,2,4); plot(y1(1:8:lim),y2(1:8:lim),'.');
axis([-8 8 -8 8]); xlabel('y1'); ylabel('y2'); title('y1-y2');

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

%---------------------------------

clear V lim z1 z2 

% end