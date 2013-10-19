function new_w = bsslearning(m1,m2,w)

l=length(m1);		% l = length(z1) = length(z2)

wm = w(1,1)*m1 + w(2,1)*m2;

g=tanh(wm);

mg1=m1'*g;
mg2=m2'*g;

Emg1=sum(mg1)/l;
Emg2=sum(mg2)/l;

gg=1./((cosh(wm)).^2);

Egg=sum(gg)/l;

new_w(1,1)=Emg1-Egg*w(1,1);
new_w(2,1)=Emg2-Egg*w(2,1);

clear Egg Ezg1 Ezg2 g gg l wz zg1 zg2