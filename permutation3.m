function  per = permutation(U1,U2,M);

% envelope����뤿����ȷ��βû��ϰ�M
% [ts-M,ts+M]���ϰϤǷ׻�,��(2M+1)��
[BIN,T]=size(U1);	% BIN:NFFT/2+1, T:����Ĺ

%-------------------------------------------
%-------------------------------------------
% Scaling���褷��ʬΥ����U

h=hilb(M+1);
h=[h(1,end:-1:2) h(1,:)];	% (2*M+1)��

% iteration
for k=1:BIN,

	% �ȷ��������ͤ�׻�
	abs1=[zeros(1,M) abs(U1(k,:)) zeros(1,M)];
	abs2=[zeros(1,M) abs(U2(k,:)) zeros(1,M)];

	for t=1:T,
	env1(t)=sum(abs1(t:t+2*M).*h);	% envelope
	env2(t)=sum(abs2(t:t+2*M).*h);
	end

	% �Ƽ��ȿ����Ȥ˼�Ǽ
	ENV1(k,:)=env1;
	ENV2(k,:)=env2;

	% ����٥��״֤���ؤ�׻�
	sim(k,1)=sum(env1.*env2)/(sqrt(sum(env1.^2))*sqrt(sum(env2.^2)));

	clear abs1 abs2 ts env1 env2	% initialize
	
end

%------------------------------------
%------------------------------------
% �����ޤǤ����Ǽ��ȿ�����¤������٥��׿������ENV������줿��
% ���줫��ºݤ�permutation��������褹�롣

% ��ؤ��㤤��˼��ȿ����¤٤�
[y,i]=sort(sim);	% y�˾��������͡����μ��ȿ�bin�ʥ�С���i

% per��permutation���Ǽ
% [0:�֤�����ɬ��̵��]
% [1:�֤�����ɬ�פ���]
per=zeros(size(sim));

% ����ͤȤ���m=1���ܤμ��ȿ�perm(m=1)�����ꡣ
m=1;		% m�ϸ��ߤη׻����
k=i(m);		% k�ϼºݤμ��ȿ�bin�ʥ�С�
per(k)=0;	% �����ܤ˷׻���Ԥ����ȿ����֤�����ɬ�פʤ����������Ȥ���

yenv1=ENV1(k,:);
yenv2=ENV2(k,:);

% i���(��ؤξ��������ȿ�bin���)permutation���褷�Ƥ���
for m=1:BIN,
	
	k=i(m);
	
	env1=ENV1(k,:);
	env2=ENV2(k,:);
	
	[a b] = max([yenv1*env1'; yenv1*env2'; yenv2*env2'; yenv2*env1']);

	if mod(b,2)==1, 		% �֤�����ɬ�פʤ�
		per(k)=0;
		yenv1=yenv1+ENV1(k,:);
		yenv2=yenv2+ENV2(k,:);
	else 	per(k)=1;		% �֤�����ɬ�פ���
		yenv1=yenv1+ENV2(k,:);
		yenv2=yenv2+ENV1(k,:);
	end

	clear a b k env1 env2 

end

% per reterned.
