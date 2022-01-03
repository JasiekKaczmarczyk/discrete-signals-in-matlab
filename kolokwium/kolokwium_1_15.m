%% zad 1
close all; clc; clear;
Fs=150;
t=0:1/Fs:6;

x=(-1*(t>=0 & t<=2)) + (-1*(t>=4 & t<=6)) + (-(3*(1-abs(t-3))).*(abs(t-3)<1));

FS=-7/3;

for n=1:80
  an=
  bn=
  
  FS=FS+an*cos(n*pi*t/3)+bn*sin(n*pi*t/3);
  
  if n==10
      x10=FS;
  end
  
  if n==80
     x80=FS; 
  end
end

plot(t,x,'b.',t,x10,'r',t,x80,'g');

%% zad 2
close all; clc;clear;

a=load('kolos/cor07.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-220:1/Fs:220;

%trojkat
tp1=-6:1/Fs:6;
x1=(1-abs(tp1)/6);
xc1=xcorr(1-x, 1-x1)+xcorr(x,x1);
nr1=find(xc1==max(xc1(:)), 1, 'first');

% piła
tp2=0:1/Fs:10;
x2=mod(1-tp2,10)/10;
xc2=xcorr(x.^5,x2.^5);
nr2=find(xc2>0.9*max(xc2(:)), 1, 'first');

% prostokat
tp3=0:1/Fs:20;
x3=0.65*ones(size(tp3));
xc3=xcorr(x,x3);
nr3=find(xc3==max(xc3(:)), 1, 'first');

%wypisanie znalezisk do kosoli
tc(nr1)
tc(nr2)
tc(nr3)

subplot(211), plot(t,x,tp3,x3);
subplot(212), plot(tc,xc3);

%%
close all; clc;clear;
Fs=400;
t=0:1/Fs:5;

xa=((1/5)*t+1).*sin(2*pi*t*45);
xb=sqrt(2)*sin(2*pi*cumsum(-2*t+20)/Fs);
xc=3*(1-abs(t-2)/0.25).*(abs(t-2)<0.25);
xd=-0.2+0.5*rand(size(t));
xe=sin(2*pi*t*5).*(t>=1 & t<=3);

x=xa+xb+xc+xd+xe;

FT=fftshift(fft(x));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

BStop=1.0*(abs(f)>=20 | abs(f)<=10);
FT_n=BStop.*FT;
x_n=ifft(ifftshift(FT_n));

figure;
subplot(211),plot(t,x,t,x_n,'r');
subplot(212),plot(f,WA,f,200*BStop,'r');

figure;
subplot(221),plot(t,xa);
subplot(222),plot(t,xb);
subplot(223),plot(t,xc);
subplot(224),plot(t,xd);