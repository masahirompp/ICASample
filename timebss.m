function B = timebss(x1,x2,r);

% xは分離対象の観測信号
% rは同時対角化する時間差の数
%----------------------------------
%----------------------------------
% whitening

% zero mean
x1=x1-mean(x1);
x2=x2-mean(x2);

% eigenvalue decomposition 'Cov'
Cov=cov(x1,x2).';	
V=eigenvalue_decomposition(Cov);

% decorrelation
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

clear r1 r2 c1 c2 m1 m2 Cov

%--------------------------------------
%--------------------------------------
% unmixing

l=length(z1);
G=[];

for tao=1:r,     %時間差tao

	m11=z1(tao+1:l).'*conj(z1(1:l-tao))/l;	%各時間差相関を計算
	m12=z1(tao+1:l).'*conj(z2(1:l-tao))/l;
	m21=z2(tao+1:l).'*conj(z1(1:l-tao))/l;
	m22=z2(tao+1:l).'*conj(z2(1:l-tao))/l;

	g=[m11-m22  m12+m21  i*(m21-m12)];
	G=[G;g];

end

[E,D]=eig(real(G'*G));
v=E(:,3);

%cov(z1,z2)
%v'*real(G'*G)*v

theta=acos(v(1,1));
cosphi=-v(2,1)/sin(theta);
sinphi=-v(3,1)/sin(theta);

G=[cos(theta/2) (cosphi+i*sinphi)*sin(theta/2);(-cosphi+i*sinphi)*sin(theta/2) cos(theta/2)];

W=G';

y1=W(1,1)*z1+W(1,2)*z2;
y2=W(2,1)*z1+W(2,2)*z2;

B=W*V;
