%%
Fs=90;
t=-2:1/Fs:2;

x=(t+1).*(t>=-2 & t<=0) + (-2*t+1).*(t>=0 & t<=1) + (-1)*(t>=1 & t<=2);

a0=-1/2;
FX=a0*ones(size(t))/2;

for n=1:80
   an=1/2*(-(2*pi*n*sin(pi*n)+4*cos(pi*n)-4)/(pi^2*n^2)-(2*pi*n*sin((pi*n)/2)+8*cos((pi*n)/2)-8)/(pi^2*n^2)-(2*(sin(pi*n)-sin((pi*n)/2)))/(pi*n));
   bn=1/2*((2*(2*sin(pi*n)-pi*n*cos(pi*n)-pi*n))/(pi^2*n^2)-(2*(4*sin((pi*n)/2)-pi*n*cos((pi*n)/2)-pi*n))/(pi^2*n^2)+(2*(cos(pi*n)-cos((pi*n)/2)))/(pi*n));
   
   FX=FX+an*cos(n*pi*t/2)+bn*sin(n*pi*t/2);
   
   if n==15
      FX15=FX; 
   end
end

plot(t,x,'.r',t,FX15,'b',t,FX,'y');

%%
close all;clc;clear;
a=load('korel_203.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-115:1/Fs:115;

tp1=-5:1/Fs:5;
x1=0.75*(1-abs(tp1)/5);
xc1=xcorr(x,x1)+xcorr(1-x,1-x1);
nr1=find(xc1==max(xc1(:)), 1, 'first');

tc(nr1)

tp2=-5:1/Fs:5;
x2=0.45*ones(size(tp2));
xc2=xcorr(x,x2)+xcorr(1-x.^5,1-x2.^5);
nr2=find(xc2==max(xc2(:)), 1, 'first');

%wypisanie czasu największej korelacji
tc(nr2)


tp3=-5:1/Fs:5;
x3=0.65*exp(-(tp3).^2/(2*2*2));
xc3=xcorr(x,x3)+xcorr(1-x.^5,1-x3.^5);
nr3=find(xc3==max(xc3(:)), 1, 'first');

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
    %N2=floor(N/2);
    %LP=exp(-(-N2:N2).^2/(2(N2/4).^2));
    %LP=LP/sum(LP(:));
    odch=k/7;
    LP=exp(-(-k:k).^2/(2*odch*odch));
    LP=LP/sum(LP(:));
    x2=conv(y,LP,'same');
    %mediana
    x3=medfilt1(y,N);
    %wienier
    x4=wiener2(y,[1,N]);

    blad(k,1)=ocena(x,x1);
    blad(k,2)=ocena(x,x2);
    blad(k,3)=ocena(x,x3);
    blad(k,4)=ocena(x,x4);
end

plot(blad);
legend;

x3=medfilt1(y,17);
ocena(x,x3)