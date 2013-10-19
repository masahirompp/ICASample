function  per = permutation(U1a,U2a,U1b,U2b,ilist);

% permutation�β��

% ɬ�פ�ʪ�����
fs=ilist(1);
NFFT=ilist(2);
shiftT=ilist(3);
T=size(U1a,2);		% ����Ĺ

% envelope����뤿����ȷ��βû��ϰ�
M=20;				% �����Ǥ�14�˸���ʹ�®���Τ����

%-------------------------------------------
%-------------------------------------------
% Scaling���褷��ʬΥ����U

% initialize
ENV1=zeros(NFFT/2+1,T);		% �Ƽ��ȿ����ȤΥ���٥��׿����׻�����ENV�˼�Ǽ����
ENV2=zeros(NFFT/2+1,T);
sim=zeros(NFFT/2+1,1);		% �Ƽ��ȿ��ˤ�����env1��env2����ؤ��Ǽ

% iteration
for k=1:NFFT/2+1,

	% �ȷ��������ͤ�׻�
	abs1=abs(U1a(k,:))+abs(U1b(k,:));
	abs2=abs(U2a(k,:))+abs(U2b(k,:));
	
	% initialize
	envelope1=zeros(1,T);
	envelope2=zeros(1,T);
	
	for ts=1:T,		

		% �ȷ��βû��ϰ�M�ˤ���ȷ��û���������
		% ts���顢�ޤ��Ͽ���Ĺ��ۤ����ϰϤ�̵��ʤ����Ǥ�0���֤��Ƥ����
		tss=[ts-M:1:ts+M];
		if ts<M, tss(1:M-ts+1)=0;
		end
		if ts>T-M,
		tss=fliplr(tss);
		tss(1:M+ts-T)=0;
		tss=fliplr(tss);
		end

		% initialize
		env1=0;
		env2=0;

		% �����ˤ������ȷ��βû�
		for m=1:length(tss)

			if tss(m)==0,env1 = env1;
			else env1 = env1 + abs1(tss(m));
			end
		
			if tss(m)==0,env2 = env2;
			else env2 = env2 + abs2(tss(m));
			end
		
		end

		% ����٥��פδ���	
		envelope1(1,ts)=env1;
		envelope2(1,ts)=env2;

		clear m env1 env2 tss
		
	end

	% �Ƽ��ȿ����Ȥ˼�Ǽ
	ENV1(k,:)=envelope1/(2*M);
	ENV2(k,:)=envelope2/(2*M);


	% ����٥��״֤���ؤ�׻�
	sim(k,1)=envelope1*envelope2'/(sqrt(envelope1*envelope1')*sqrt(envelope2*envelope2'));

	clear abs1 abs2 ts envelope1 envelope2
	
end

%------------------------------------
%------------------------------------
% �����ޤǤ����Ǽ��ȿ�����¤������٥��׿������ENV������줿��
% ���줫��ºݤ�permutation��������褹�롣

% �ޤ���permutation����ꤷ�Ƥ����ȿ�bin�ν��֤���롣

minbin=zeros(size(sim));	% ��ؤξ��������ȿ��ʥ�С����˼�Ǽ

for m=1:NFFT/2+1,

	[a b]=min(sim);

	minbin(m)=b;
	sim(b)=100;

	clear a b

end
%minbin'
% minbin�ν��֤���ؤ�׻�����permutation���衣

% per��permutation���Ǽ
% [0:�֤�����ɬ��̵��]
% [1:�֤�����ɬ�פ���]
per=zeros(size(sim));

% ����ͤȤ���m=1���ܤμ��ȿ�perm(m=1)�����ꡣ
m=1;		% m�ϸ��ߤη׻����
k=minbin(m);	% k�ϼºݤμ��ȿ�bin�ʥ�С�
per(k)=0;	% �����ܤ˷׻���Ԥ����ȿ����֤�����ɬ�פʤ����������Ȥ���

yenv1=ENV1(k,:);
yenv2=ENV2(k,:);

% m=2���ܰʹߤμ��ȿ�������ؤǷ��ꤷ�Ƥ���
for m=2:NFFT/2+1,
	
	k=minbin(m);
	
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
