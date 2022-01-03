%CZAS 24:42 min

%%
close all; clc; clear;

Fs=120;
t=0:1/Fs:6;

x=1*((t>=0 & t<=2) | (t>=4 & t<=6)) + 4*(1-abs(t-3)).*(abs(t-3)<1);

a0=8/3;
FX=a0*ones(size(t))/2;

for n=1:80
    an=1/3*(-(36*(cos((4*pi*n)/3)-2*cos(pi*n)+cos((2*pi*n)/3)))/(pi^2*n^2)+(3*sin((2*pi*n)/3))/(pi*n)+(3*(sin(2*pi*n)-sin((4*pi*n)/3)))/(pi*n));
    bn=1/3*(-(36*(sin((4*pi*n)/3)-2*sin(pi*n)+sin((2*pi*n)/3)))/(pi^2*n^2)-(3*cos((2*pi*n)/3)-3)/(pi*n)-(3*(cos(2*pi*n)-cos((4*pi*n)/3)))/(pi*n));
   
   FX=FX+an*cos(n*pi*t/3)+bn*sin(n*pi*t/3);
   
   if n==15
       FX15=FX;
   end
   
end

plot(t,x,'.b',t,FX15,'r',t,FX,'g');

%%
close all;clc;clear;
a=load('kolos/cor06.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-55:1/Fs:55;

tp1=-1.5:1/Fs:1.5;
x1=(1-abs(tp1)/1.5);
xc1=xcorr(1-x,1-x1)+xcorr(x,x1);
nr1=find(xc1==max(xc1(:)), 1, 'first');

tp2=-3:1/Fs:3;
x2=0.9*exp(-(tp2).^2/(2*0.75^2));
xc2=xcorr(sqrt(1-x),sqrt(1-x2))+xcorr(sqrt(x),sqrt(x2));
nr2=find(xc2==max(xc2(:)), 3, 'first');

tp3=0:1/Fs:5;
x3=0.65*ones(size(tp3));
xc3=xcorr(x,x3)+xcorr(1-x.^5,1-x3.^5);
nr3=find(xc3==max(xc3(:)), 1, 'first');

tc(nr1)
tc(nr2)
tc(nr3)

subplot(211),plot(t,x,tp1,x1,'r',tp2,x2,'k',tp3,x3,'g');
subplot(212),plot(tc,xc3);

%%

close all;clc;clear;
Fs=250;
t=0:1/Fs:8;

xa=2*sqrt(2)*exp(-(t-4.5).^2/(2*0.15^2));
xb=2.5*sin(2*pi*t.*(1/4*t+11));
xc=1.1*sin(2*pi*t*7.5).*(t>=1 & t<=4);
xd=(1/16*t+1.7).*sin(2*pi*t*39);
xe=0.15*randn(size(t));

x=xa+xb+xc+xd+xe;

FT=fftshift(fft(x));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

BStop=1.0*(abs(f)>=40 | abs(f)<38.5);
FT_n=BStop.*FT;
x_n=ifft(ifftshift(FT_n));

figure;
subplot(211),plot(t,x,t,x_n,'r');
subplot(212),plot(f,WA,f,600*BStop,'r');

figure;
subplot(221),plot(t,xa);
subplot(222),plot(f,xb);
subplot(223),plot(t,xc);
subplot(224),plot(f,xd);