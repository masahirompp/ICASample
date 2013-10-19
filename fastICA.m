function new_w = fastICA(w,z1,z2)

%fast ICA algorithm

y=(conj(w(1))).*z1 + (conj(w(2))).*z2;
yy=y.*conj(y);					% yyは実数
gyy=((yy+0.1).^(-1/2))./2;			% 実数
ggyy=-((yy+0.1).^(-3/2))./4;			% 実数

E1=mean(z1.*conj(y).*gyy);
E2=mean(z2.*conj(y).*gyy);

E3=mean(gyy+yy.*ggyy);

new_w=[E1-E3*w(1); E2-E3*w(2)];

clear y yy gyy ggyy E1 E2 E3
