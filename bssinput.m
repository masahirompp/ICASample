% Script file : bssinput.m

% return x1, x2, fs, NFFT, shiftT, learntime, ilist

%-------------------

x1=input('Enter a observed signal x1: ');
x2=input('Enter the other observed signal x2: ');

sa=length(x1)-length(x2);
if sa>0, x2 = [x2; zeros(sa,1)];
elseif sa<0, x1=[x1; zeros(-sa,1)];                     
end

% prevent culcuration error
x1=10*x1;
x2=10*x2;

%-----------------

ilist=input('Enter a input list vector, ilist, [fs; NFFT; shiftT; learntime(sec)]: ');

if length(ilist)~=4
	disp('lack of object');
	ilist=input('Enter a input list vector, ilist, [fs; NFFT; shiftT; learntime(sec)]: ');
end

fs=ilist(1);
NFFT=ilist(2);
shiftT=ilist(3);
learntime=ilist(4);

%------------------

if log2(NFFT) ~= round(log2(NFFT)),
	disp('NFFT = 2^k.');
	learntime=input('N-point FFT, NFFT: ');
end

if learntime > length(x1)/fs,
	disp('The learntime is longer than the length of x1 and x2.');
	learntime=input('learning and unmixing time [sec] : ');
end

clear sa

