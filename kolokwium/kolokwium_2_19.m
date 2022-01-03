%%
Fs=125;
t=-3:1/Fs:1;

x=(t+1).*(t>=-3 & t<=-1) + (2*t+2).*(t>=-1 & t<=0) + (1)*(t>=0 & t<=1);

a0=0;
FX=zeros(size(t));

for n=1:75
   an=1/2*(-(4*(pi*n*sin((3*pi*n)/2)+cos((3*pi*n)/2)-cos((pi*n)/2)))/(pi^2*n^2)-(8*cos((pi*n)/2)-8)/(pi^2*n^2)+(2*sin((pi*n)/2))/(pi*n));
   bn=1/2*((4*(sin((3*pi*n)/2)-pi*n*cos((3*pi*n)/2)-sin((pi*n)/2)))/(pi^2*n^2)+(8*sin((pi*n)/2))/(pi^2*n^2)-4/(pi*n)-(2*cos((pi*n)/2)-2)/(pi*n));
   
   FX=FX+an*cos(n*pi*t/2)+bn*sin(n*pi*t/2);
   
   if n==18
      FX18=FX; 
   end
end

plot(t,x,'.b',t,FX18,'g',t,FX,'k');

%%
close all;clc;clear;
a=load('korel_206.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-222:1/Fs:222;

tp1=-10:1/Fs:10;
x1=1.7*(1-abs(tp1)/10);
xc1=xcorr(x,x1)+xcorr(1-x,1-x1);
nr1=find(xc1==max(xc1(:)), 2, 'first');

tc(nr1)

tp2=-7:1/Fs:7;
x2=1.2*ones(size(tp2));
xc2=xcorr(1-x.^3,1-x2.^3)+xcorr(x,x2);
nr2=find(xc2==max(xc2(:)), 1, 'first');

%wypisanie czasu największej korelacji
tc(nr2)


tp3=-20:1/Fs:20;
x3=1.2*exp(-(tp3).^2/(2*5^2));
xc3=xcorr(x.^2,x3.^2);
nr3=find(xc3>min(xc3(:))+0.2, 1, 'first');

tc(nr3)

subplot(211),plot(t,x,tp1,x1,'g',tp2,x2,'r',tp3,x3,'y');
subplot(212),plot(tc,xc3);

%%
close all;clc;clear;
a=load('cad_201.txt');
x=a(:,1)';
y=a(:,2)';
Fs=500;
t=linspace(0,length(x)/Fs,length(x));

ocena=@(x,y) sqrt((1/length(t)*sum((x(:)-y(:)).^2)));

FTx=fftshift(fft(x));
WAx=abs(FTx);

FTy=fftshift(fft(y));
WAy=abs(FTy);
f=linspace(-Fs/2,Fs/2,length(t));

f0=15;
n=4;
Butterworth=(1.0)./(1+(f/f0).^(2*n)).*(abs(f)<=22);
Y=real(ifft(ifftshift(Butterworth.*FTy)));

L2=ocena(x,Y);

subplot(211),plot(t,x,'b',t,Y,'r');
subplot(212), plot(f,WAx,'r',f,WAy,'b',f,10000*Butterworth,'k');
L2

%%
close all;clc;clear;
a=load('cad_201.txt');
x=a(:,1)';
y=a(:,2)';
Fs=500;
t=linspace(0,length(x)/Fs,length(x));

ocena=@(x,y) sqrt((1/length(t)*sum((x(:)-y(:)).^2)));
blad=zeros(50,4);

for k=1:50
    N=2*k+1;
    %zwykly filtr
    x1=conv(y,ones(1,N)/N, 'same');
    % gauss
    N2=floor(N/2);
    LP=exp(-(-N2:N2).^2/(2*(N2/4).^2));
    LP=LP/sum(LP(:));
    x2=conv(y,LP,'same');
    %mediana
    x3=medfilt1(y,N);
    %wienier
    x4=wiener2(y,[1,N]);
    
    x5=wiener2(medfilt1(wiener2(y,[1,N]),N),[1,N]);
    x6=wiener2(wiener2(medfilt1(y,N),[1,N]),[1,N]);
    x7=medfilt1(medfilt1(medfilt1(y,N),N),N);

    blad(k,1)=ocena(x,x3);
    blad(k,2)=ocena(x,x4);
    blad(k,3)=ocena(x,x5);
    blad(k,4)=ocena(x,x7);
end

%plot(blad);
%legend;

FTy=fftshift(fft(y));
WAy=abs(FTy);

FTx=fftshift(fft(x));
WAx=abs(FTx);

f=linspace(-Fs/2,Fs/2,length(t));

f0=4;
n=1;
Butterworth=(1.0)./(1+(f/f0).^(2*n));
c=real(ifft(ifftshift(Butterworth.*FTy)));

N=2*4+1;
xdobry=wiener2(medfilt1(wiener2(c,[1,N]),N),[1,N]);

ocena(x,xdobry)

subplot(211),plot(t,x,'b');
subplot(212), plot(f,WAx,'b',f,WAy,'k',f,10000*Butterworth,'r');