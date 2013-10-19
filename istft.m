function [y]=istft(Xk,nfft,fs,shift_size)
%   Inverse Short Time Fourier Transfrom
%                     date: 29 Nov 2007, Shinpei Tsuchiya
%  <Usage>
%  y=istft(Xk,nfft,fs,shift_size)
%  Xk : stft matrix(nfft x number of frame)
%  y  : synthesis time signal (Longitudinal Vector)

number_of_frame=size(Xk,2);
length_of_signal=shift_size*(number_of_frame-1)+nfft;
y=zeros(length_of_signal,1);
xn=real(ifft(Xk,nfft));
for frame=1:number_of_frame
       offset=shift_size*(frame-1);
       y(offset+1:offset+nfft)=y(offset+1:offset+nfft)+xn(:,frame).*hanning(nfft);
end



