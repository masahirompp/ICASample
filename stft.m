function Xk=stft(y,nfft,fs,shift_size)
%   Short Time Fourier Transfrom
%                  date: 29 Nov 2007, Shinpei Tsuchiya
%  <Usage>
%  [Xk]=stft(y,nfft,fs,shift_size)
%  y  : time signal (Longitudinal Vector)
%  Xk : stft matrix (nfft x number_of_frame)
% +++ default window function is hanning window +++


length_of_signal=length(y);
number_of_frame=floor((length_of_signal-(nfft-shift_size))/shift_size);
Xk=zeros(nfft,number_of_frame);
F=[0:nfft-1]*fs/nfft;
T=shift_size/fs/2*[1:number_of_frame];
window_f=hanning(nfft); % default window funciton is hanning window

for frame=1:number_of_frame
    offset=shift_size*(frame-1);
    xn=hanning(nfft).*y(offset+1:offset+nfft);
    Xk(:,frame)=fft(xn,nfft);
end

 %imagesc(T,F,20*log10(abs(Xk)));
 %axis xy;grid,ylim([0 fs/2]);
 %xlabel('time(s)');
 %ylabel('frequency(Hz)');
