%% zad3
close all;clc;clear;
x=load('kolos/szum07.txt');
a=x(:,1)';
b=x(:,2)';
Fs=200;
t=linspace(0,length(a)/Fs,length(a));

ocena=@(x,y) (1/length(t)*sum(abs(x(:)-y(:))));
blad=zeros(49,4);

FTa=fftshift(fft(a));
WAa=abs(FTa);

FTb=fftshift(fft(b));
WAb=abs(FTb);

f=linspace(-Fs/2,Fs/2,length(t));

f0=3.6;
n=1;
Butterworth=(1.0)./(1+(f/f0).^(2*n)).*(abs(f)>=0 & abs(f)<=13);
c=real(ifft(ifftshift(Butterworth.*FTb)));

ocena(a,c)

subplot(211),plot(t,a,'b');
subplot(212), plot(f,WAa,f,WAb,'r',f,100000*Butterworth,'g');