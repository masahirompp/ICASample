

tx=(1:length(x1))/fs;
ty=(1:length(y1))/fs;
figure;
subplot(2,3,1);plot(tx,ch1);axis([0,length(y1)/fs,-1,1]);title('ch1')
subplot(2,3,2);plot(tx,x1);axis([0,length(y1)/fs,-1,1]);title('x1')
subplot(2,3,3);plot(ty,y1);axis([0,length(y1)/fs,-1,1]);title('y1')
subplot(2,3,4);plot(tx,ch2);axis([0,length(y1)/fs,-1,1]);title('ch2')
subplot(2,3,5);plot(tx,x2);axis([0,length(y1)/fs,-1,1]);title('x2')
subplot(2,3,6);plot(ty,y2);axis([0,length(y1)/fs,-1,1]);title('y2')

figure;
subplot(2,3,1);specgram(ch1,1024,16000,hann(1024),1000);title('ch1')
subplot(2,3,2);specgram(x1,1024,16000,hann(1024),1000);title('x1')
subplot(2,3,3);specgram(y1,1024,16000,hann(1024),1000);title('y1')
subplot(2,3,4);specgram(ch2,1024,16000,hann(1024),1000);title('ch2')
subplot(2,3,5);specgram(x2,1024,16000,hann(1024),1000);title('x2')
subplot(2,3,6);specgram(y2,1024,16000,hann(1024),1000);title('y2')

a1=ch1(40001:41024);
a1=a1-mean(a1);
a2=ch2(40001:41024);
a2=a2-mean(a2);
b1=x1(40001:41024);
b1=b1-mean(b1);
b2=x2(40001:41024);
b2=b2-mean(b2);
c1=y1(40001:41024);
c1=c1-mean(c1);
c2=y2(40001:41024);
c2=c2-mean(c2);

figure;
subplot(2,1,1);hold on;psd(a1,1024,16000,hann(1024))
psd(b1,1024,16000,hann(1024));hold off;title('filtering ch1 x1')
subplot(2,1,2);hold on;psd(a2,1024,16000,hann(1024))
psd(b2,1024,16000,hann(1024));hold off;title('filtering ch2 x2')

figure;
subplot(2,1,1);hold on;psd(x1,1024,16000,hann(1024))
psd(y1,1024,16000,hann(1024));
psd(y2,1024,16000,hann(1024));hold off
title('bss x1 y1 y2')
subplot(2,1,2);hold on;psd(x2,1024,16000,hann(1024))
psd(y1,1024,16000,hann(1024));
psd(y2,1024,16000,hann(1024));hold off
title('bss x1 y1 y2')
