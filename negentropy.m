function J = negentropy(y1,y2,l,Hu)

% G and H is non-liner function

G1=log(cosh(y1));
G2=log(cosh(y2));
H1=-exp((-y1.^2)/2);
H2=-exp((-y2.^2)/2);

J1= ((sum(G1))/l)^2 + ((sum(H1))/l - (sum(Hu))/l)^2;
J2= ((sum(G2))/l)^2 + ((sum(H2))/l - (sum(Hu))/l)^2;

J=max(J1,J2);

clear y1 y2 G1 G2 H1 H2 