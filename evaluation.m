function E1 = evaluation(A)

% text p308
%
% P=[p11 p12; p21 p22]
% m1=[p11 p12]
% m2=[p21 p22]
% n1=[p11 p21]
% n2=[p12 p22]

P=abs(A);

m1=P(1,:);
m2=P(2,:);
n1=P(:,1)';
n2=P(:,2)';

E1=sum(m1)/max(m1)+sum(m2)/max(m2)+sum(n1)/max(n1)+sum(n2)/max(n2)-4;

clear m1 m2 n1 n2 A