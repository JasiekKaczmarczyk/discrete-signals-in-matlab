%%
close all;clc;clear;
a=load('dem_203.txt');
x=a(:,1)';
W=a(:,2)';
Fs=500;
t=linspace(0,length(x)/Fs,length(x));

ocena=@(x,y) sqrt((1/length(t)*sum((x(:)-y(:)).^2)));
blad=zeros(50,4);

for k=1:50
    N=2*k+1;
    %zwykly filtr
    x1=conv(W,ones(1,N)/N, 'same');
    % gauss
    N2=floor(N/2);
    LP=exp(-(-N2:N2).^2/(2*(N2/4).^2));
    LP=LP/sum(LP(:));
    x2=conv(W,LP,'same');
    %mediana
    x3=medfilt1(W,N);
    %wienier
    x4=wiener2(W,[1,N]);
    
    x5=wiener2(medfilt1(wiener2(W,[1,N]),N),[1,N]);
    x6=wiener2(wiener2(medfilt1(W,N),[1,N]),[1,N]);

    blad(k,1)=ocena(x,x1);
    blad(k,2)=ocena(x,x2);
    blad(k,3)=ocena(x,x5);
    blad(k,4)=ocena(x,x6);
end

%k=15 0.55
N=2*7+1;
xdobry=wiener2(medfilt1(wiener2(W,[1,N]),N),[1,N]);
ocena(x,xdobry)

%
FTx=fftshift(fft(x));
WAx=abs(FTx);

FTw=fftshift(fft(W));
WAw=abs(FTw);
f=linspace(-Fs/2,Fs/2,length(t));

f0=1;
n=1;
Butterworth=(1.0)./(1+(f/f0).^(2*n));
Y=real(ifft(ifftshift(Butterworth.*FTw)));

ocena(x,Y)

subplot(211),plot(t,x,'b',t,Y,'r');
subplot(212), plot(f,WAx,'r',f,WAw,'b',f,10000*Butterworth,'k');
%plot(blad);
%legend;