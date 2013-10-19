function Z1= test(X1,NFFT)

ans=size(X1);
Z1=zeros(1,ans(2));

bin=1;

while bin<NFFT/2+1

	bin=bin+1;

	x1=X1(bin,:);

	x1=x1./((var(x1))^(1/2));
	
	Z1=[Z1;x1];

end

