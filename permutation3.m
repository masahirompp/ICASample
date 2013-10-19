function  per = permutation(U1,U2,M);

% envelopeを求めるための波形の加算範囲M
% [ts-M,ts+M]の範囲で計算,全(2M+1)点
[BIN,T]=size(U1);	% BIN:NFFT/2+1, T:信号長

%-------------------------------------------
%-------------------------------------------
% Scalingを解決した分離信号U

h=hilb(M+1);
h=[h(1,end:-1:2) h(1,:)];	% (2*M+1)点

% iteration
for k=1:BIN,

	% 波形の絶対値を計算
	abs1=[zeros(1,M) abs(U1(k,:)) zeros(1,M)];
	abs2=[zeros(1,M) abs(U2(k,:)) zeros(1,M)];

	for t=1:T,
	env1(t)=sum(abs1(t:t+2*M).*h);	% envelope
	env2(t)=sum(abs2(t:t+2*M).*h);
	end

	% 各周波数ごとに収納
	ENV1(k,:)=env1;
	ENV2(k,:)=env2;

	% エンベロープ間の相関を計算
	sim(k,1)=sum(env1.*env2)/(sqrt(sum(env1.^2))*sqrt(sum(env2.^2)));

	clear abs1 abs2 ts env1 env2	% initialize
	
end

%------------------------------------
%------------------------------------
% ここまでの操作で周波数順に並んだエンベロープ信号行列ENVが得られた。
% これから実際にpermutationの問題を解決する。

% 相関の低い順に周波数を並べる
[y,i]=sort(sim);	% yに昇順に相関値、その周波数binナンバーがi

% perにpermutationを収納
% [0:置き換え必要無し]
% [1:置き換え必要あり]
per=zeros(size(sim));

% 初期値としてm=1番目の周波数perm(m=1)を設定。
m=1;		% mは現在の計算回数
k=i(m);		% kは実際の周波数binナンバー
per(k)=0;	% 一番目に計算を行う周波数は置き換え必要なし、これを基準とする

yenv1=ENV1(k,:);
yenv2=ENV2(k,:);

% i順に(相関の小さい周波数bin順に)permutationを解決していく
for m=1:BIN,
	
	k=i(m);
	
	env1=ENV1(k,:);
	env2=ENV2(k,:);
	
	[a b] = max([yenv1*env1'; yenv1*env2'; yenv2*env2'; yenv2*env1']);

	if mod(b,2)==1, 		% 置き換え必要なし
		per(k)=0;
		yenv1=yenv1+ENV1(k,:);
		yenv2=yenv2+ENV2(k,:);
	else 	per(k)=1;		% 置き換え必要あり
		yenv1=yenv1+ENV2(k,:);
		yenv2=yenv2+ENV1(k,:);
	end

	clear a b k env1 env2 

end

% per reterned.
