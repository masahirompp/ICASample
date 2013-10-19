% Script file : bssmix.m

s1=input('Enter a original signal s1: ');
s2=input('Enter the other original signal s2: ');

% make a same length signal
sa=length(s1)-length(s2);
if sa>0, s2 = [s2; zeros(sa,1)];
elseif sa<0, s1=[s1; zeros(-sa,1)];                     
end
clear sa;

% make mixing matrix
A=tanh(2*randn(2));

% make observed singnal x
x1=A(1,1)*s1+A(1,2)*s2;
x2=A(2,1)*s1+A(2,2)*s2;

disp(' ');
disp('mixing matrix A');
A
