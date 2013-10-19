function TEST = test2(X1,X2,NFFT)

X1=test(X1,NFFT);
X2=test(X2,NFFT);

bin=1;
TEST=[];

while bin < NFFT/2

	bin=bin+1;

	m2 = X1(bin,:);
	m3 = X1(bin+1,:);

	n2 = X2(bin,:);
	n3 = X2(bin+1,:);

	l1 = abs(m2*m3'/length(m2));
	l2 = abs(m2*n3'/length(m2));
	l3 = abs(n2*m3'/length(m2));
	l4 = abs(n2*n3'/length(m2));

	if (l1>l2)+(l3<l4)==0, k=2;
	elseif (l1>l2)+(l3<l4)==1, k=1;
	elseif (l1>l2)+(l3<l4)==2, k=0;
	else k=3;
	end

	ans=[l1 l2 l3 l4];
	TEST=[TEST; bin ans k];

end

