function  per = permutation(U1a,U2a,U1b,U2b,ilist);

% permutationの解決

% 必要な物を準備
fs=ilist(1);
NFFT=ilist(2);
shiftT=ilist(3);
T=size(U1a,2);		% 信号長

% envelopeを求めるための波形の加算範囲
M=20;				% ここでは14に固定（高速化のため）

%-------------------------------------------
%-------------------------------------------
% Scalingを解決した分離信号U

% initialize
ENV1=zeros(NFFT/2+1,T);		% 各周波数ごとのエンベロープ信号を計算し、ENVに収納する
ENV2=zeros(NFFT/2+1,T);
sim=zeros(NFFT/2+1,1);		% 各周波数におけるenv1とenv2の相関を収納

% iteration
for k=1:NFFT/2+1,

	% 波形の絶対値を計算
	abs1=abs(U1a(k,:))+abs(U1b(k,:));
	abs2=abs(U2a(k,:))+abs(U2b(k,:));
	
	% initialize
	envelope1=zeros(1,T);
	envelope2=zeros(1,T);
	
	for ts=1:T,		

		% 波形の加算範囲Mによる波形加算の前設定
		% tsが負、または信号長を越える範囲は無視（ここでは0と置いている）
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

		% 各点における波形の加算
		for m=1:length(tss)

			if tss(m)==0,env1 = env1;
			else env1 = env1 + abs1(tss(m));
			end
		
			if tss(m)==0,env2 = env2;
			else env2 = env2 + abs2(tss(m));
			end
		
		end

		% エンベロープの完成	
		envelope1(1,ts)=env1;
		envelope2(1,ts)=env2;

		clear m env1 env2 tss
		
	end

	% 各周波数ごとに収納
	ENV1(k,:)=envelope1/(2*M);
	ENV2(k,:)=envelope2/(2*M);


	% エンベロープ間の相関を計算
	sim(k,1)=envelope1*envelope2'/(sqrt(envelope1*envelope1')*sqrt(envelope2*envelope2'));

	clear abs1 abs2 ts envelope1 envelope2
	
end

%------------------------------------
%------------------------------------
% ここまでの操作で周波数順に並んだエンベロープ信号行列ENVが得られた。
% これから実際にpermutationの問題を解決する。

% まずはpermutationを決定してく周波数binの順番を決める。

minbin=zeros(size(sim));	% 相関の小さい周波数ナンバーを順に収納

for m=1:NFFT/2+1,

	[a b]=min(sim);

	minbin(m)=b;
	sim(b)=100;

	clear a b

end
%minbin'
% minbinの順番に相関を計算し、permutationを解決。

% perにpermutationを収納
% [0:置き換え必要無し]
% [1:置き換え必要あり]
per=zeros(size(sim));

% 初期値としてm=1番目の周波数perm(m=1)を設定。
m=1;		% mは現在の計算回数
k=minbin(m);	% kは実際の周波数binナンバー
per(k)=0;	% 一番目に計算を行う周波数は置き換え必要なし、これを基準とする

yenv1=ENV1(k,:);
yenv2=ENV2(k,:);

% m=2番目以降の周波数を順に相関で決定していく
for m=2:NFFT/2+1,
	
	k=minbin(m);
	
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
