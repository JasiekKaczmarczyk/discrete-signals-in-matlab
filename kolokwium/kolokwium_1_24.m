%%
close all;clc;clear;

Fs=150;
t=0:1/Fs:6;

x=-1*((t>=0 & t<=2) | (t>=4 & t<=6)) + (3*(1-abs(t-3))-2).*(abs(t-3)<1);

a0=-5/6;
FX=a0*ones(size(t));

for n=1:120
    %f parzysta
   an=1/3*(-(6*pi*n*sin((4*pi*n)/3)+27*cos((4*pi*n)/3)-54*cos(pi*n)-6*pi*n*sin((2*pi*n)/3)+27*cos((2*pi*n)/3))/(pi^2*n^2)-(3*sin((2*pi*n)/3))/(pi*n)-(3*(sin(2*pi*n)-sin((4*pi*n)/3)))/(pi*n));
   
   FX=FX+an*cos(n*pi*t/3);
   
   if n==8
      FX8=FX; 
   end
end

plot(t,x,'.b',t,FX8,'r',t,FX,'g');

%%
close all;clc;clear;
a=load('cor06.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-max(t):1/Fs:max(t);

tp1=-1.5:1/Fs:1.5;
x1=(1-abs(tp1)/1.5);
xc1=xcorr(x,x1)+xcorr(1-x,1-x1);
nr1=find(xc1==max(xc1(:)), 1, 'first');

tc(nr1)

tp2=-2.25:1/Fs:2.25;
x2=0.9*exp(-(tp2).^2/(2*0.75^2));
xc2=xcorr(sqrt(x),sqrt(x2))+xcorr(sqrt(1-x),sqrt(1-x2));
nr2=find(xc2==max(xc2(:)), 3, 'first');

tc(nr2)+2.25


tp3=-2.5:1/Fs:2.5;
x3=0.65*ones(size(tp3));
xc3=xcorr(x,x3);
nr3=find(xc3>0.9*max(xc3(:)), 1, 'first');

tc(nr3)

subplot(211),plot(t,x,tp1,x1,'g',tp2,x2,'r',tp3,x3,'y');
subplot(212),plot(tc,xc2);

%% zad3
close all;clc;clear;

Fs=400;
t=0:1/Fs:6;

xa=sqrt(3)*sin(2*pi*t.*(-t+22));
xb=(-1/6*t+3).*sin(2*pi*t*43);
xc=2*(1-abs(t-3)/0.3).*(abs(t-3)<0.3);
xd=sin(2*pi*t*6).*(t>=2 & t<=4);
xe=-0.3+(0.2-(-0.3))*rand(size(t));

x=xa+xb+xc+xd+xe;

FT=fftshift(fft(x));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

BStop=1.0*(abs(f)>=22 | abs(f)<=10);
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