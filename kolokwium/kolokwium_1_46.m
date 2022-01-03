%%
close all; clc; clear;

Fs=120;
t=0:1/Fs:6;

x=-1*((t>0 & t<2) | (t>4 & t<6)) - 4*(1-abs(t-3)).*(abs(t-3)<1);

a0=-8/3;
FX=a0*ones(size(t))/2;

for n=1:85
   an=1/3*((36*(cos((4*pi*n)/3)-2*cos(pi*n)+cos((2*pi*n)/3)))/(pi^2*n^2)-(3*sin((2*pi*n)/3))/(pi*n)-(3*(sin(2*pi*n)-sin((4*pi*n)/3)))/(pi*n));
   bn=1/3*((36*(sin((4*pi*n)/3)-2*sin(pi*n)+sin((2*pi*n)/3)))/(pi^2*n^2)+(3*cos((2*pi*n)/3)-3)/(pi*n)+(3*(cos(2*pi*n)-cos((4*pi*n)/3)))/(pi*n));
   
   FX=FX+an*cos(n*pi*t/3)+bn*sin(n*pi*t/3);
   
   if n==14
       FX14=FX;
   end
   
end

plot(t,x,'.b',t,FX14,'k',t,FX,'y');

%%
close all;clc;clear;
a=load('kolos/cor05.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-110:1/Fs:110;

%trojkat
%czas probki
tp1=-3:1/Fs:3;
%tworzenie trojkata
x1=1*(1-abs(tp1)/3);
%korelacja
xc1=xcorr(x,x1)+xcorr(1-x,1-x1);
%znajdowanie największej korelacji
nr1=find(xc1==max(xc1(:)), 1, 'first');

%wypisanie czasu największej korelacji
tc(nr1)

%gauss
%czas probki
tp2=-5:1/Fs:5;
%tworzenie trojkata
x2=0.9*exp(-(tp2).^2/(2*1.5^2));
%korelacja
xc2=xcorr(x.^5,x2.^5);
%znajdowanie największej korelacji
nr2=find(xc2==max(xc2(:)), 3, 'first');

%wypisanie czasu największej korelacji
tc(nr2)

%trojkat
%czas probki
tp3=-5:1/Fs:5;
%tworzenie trojkata
x3=0.65*ones(size(tp3));
%korelacja
xc3=xcorr(x,x3);
%znajdowanie największej korelacji
nr3=find(xc3==max(xc3(:)), 1, 'first');

tc(nr3)

subplot(211),plot(t,x,tp2,x2);
subplot(212),plot(tc,xc2);

%%
close all;clc;clear;
Fs=300;
t=0:1/Fs:8;

xa=(-0.7/8*t+2.5).*sin(2*pi*t*37);
xb=2*exp(-(t-3).^2/(2*0.2^2));
xc=sqrt(5)*sin(2*pi*cumsum(t+10)/Fs);
xd=sin(2*pi*t*7).*(t>=2 & t<=4);
xe=0.1*randn(size(t));

x=xa+xb+xc+xd+xe;

FT=fftshift(fft(x));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

BStop=1.0*(abs(f)>=18 | abs(f)<10);
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