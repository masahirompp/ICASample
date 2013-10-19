function WV = unmixing(x1,x2);

%----------------------------
%----------------------------
% pre-whitening

% zero mean
x1=x1-mean(x1);
x2=x2-mean(x2);

% decorrelation
Cov=cov(x1,x2).';
V=eigenvalue_decomposition(Cov);
z1=V(1,1)*x1+V(1,2)*x2;
z2=V(2,1)*x1+V(2,2)*x2; 

% normalize var(real(z))=0.5, var(imag(z))=0.5, var(z)=1

r1=real(z1);
r2=real(z2);
c1=imag(z1);
c2=imag(z2);

r1=r1./sqrt(2*var(r1));
r2=r2./sqrt(2*var(r2));

if var(c1)==0,c1=0;		% bin=bin_numberでは虚数部分が0であるため、分散も0になる
else c1=c1./sqrt(2*var(c1));	% そのために値が発散するのを防ぐ
end

if var(c2)==0,c2=0;
else c2=c2./sqrt(2*var(c2));
end

z1=r1+c1*i;
z2=r2+c2*i;

% z1,z2,V returned

%----------------------------
%----------------------------
% unmixing

%------------------
% initialize w1 w2

W=randn(2,2)+randn(2,2)*i;
w1=W(1,:).';
w2=W(2,:).';
w1=w1/norm(w1);
w2=w2/norm(w2);
W=[w1 w2].';

U=eigenvalue_decomposition(W*W');
W=U*W;

%------------------
% fast ICA algorithm

for counter=1:5,

	w1=W(1,:).';
	w2=W(2,:).';

	% learning
	w1=fastICA(w1,z1,z2);
	w2=fastICA(w2,z1,z2);
	w1=w1/norm(w1);
	w2=w2/norm(w2);

	% orthogonalization
	W=[w1 w2].';
	U=eigenvalue_decomposition(W*W');
	W=U*W;

end

WV=W*V;			% 分離行列B

y1=w1(1)*z1+w1(2)*z2;
y2=w2(1)*z1+w2(2)*z2;
covy=cov(y1,y2).'
