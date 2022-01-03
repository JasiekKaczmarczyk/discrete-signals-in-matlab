%% zad 1
close all; clc; clear;

Fs=200;
t=-1:1/Fs:3;

x=-(2*(1-abs(t))-1).*(abs(t)<1)+1.0*(t>=1);

a0=2/4;
FX=a0*ones(size(t));

for n=1:80
   %funkcja nie jest ani parzysta ani nieparzysta
   an=1/2*((4*pi*n*sin((pi*n)/2)+16*cos((pi*n)/2)-16)/(pi^2*n^2)+(2*(sin((3*pi*n)/2)-sin((pi*n)/2)))/(pi*n));
   bn=1/2*(-(2*(cos((3*pi*n)/2)-cos((pi*n)/2)))/(pi*n));
   
   FX=FX+an*cos(n*pi*t/2)+bn*sin(n*pi*t/2);
   
   if n==25
       FX25=FX;
   end
   
end

plot(t,x,'.k',t,FX25,'y',t,FX,'r');

%% zad 2
close all; clc; clear;

Fs=40;
t=0:1/Fs:50;

xa=sqrt(3)/2*sin(2*pi*t.*(-3.5/80*t+6));
xb=0.6*(1-abs(t-24)/10).*(abs(t-24)<10);
xc=(0.3/40*t+0.6).*sin(2*pi*t*9);
xd=0.9*mod(t,0.8)/0.8;

%szum impulsowy
x1=0.85*exp(-(t.^2)/(2*0.02^2));
%trojkat
x2=0.5*(1-abs(t-0.4)/0.2).*(abs(t-0.4)<0.2);

xe=conv(x1,x2,'same');

x=xa+xb+xc+xd+xe;

FT=fftshift(fft(x));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

LP=1.0*(abs(f)<=5);
f0=1.5;
GaussS=exp(-f.^2/(2*f0^2));
Gauss=conv(LP,GaussS,'same');
x_n=ifft(ifftshift(Gauss.*FT));

figure;
subplot(211),plot(t,x,'b',t,x_n,'r');
subplot(212),plot(f,WA,'b',f,Gauss,'r');

figure;
subplot(221),plot(t,xa);
subplot(222),plot(t,xb);
subplot(223),plot(t,xc);
subplot(224),plot(t,xd);

%% zad 3
close all;clc;clear;
a=load('szum_108.txt');
%niezaszumiony
x=a(:,1)';
%zaszumiony
xs=a(:,2)';
Fs=1000;
t=linspace(0,length(x)/Fs,length(x));

ocena=@(x,y) sqrt((1/length(t)*sum((x(:)-y(:)).^2)));
blad=zeros(30,4);

for k=1:30
    N=2*k+1;
    
    x1=conv(xs,ones(1,N)/N, 'same');
    
    N2=floor(N/2);
    LP=exp(-(-N2:N2).^2/(2*(N2/4).^2));
    LP=LP/sum(LP(:));
    x2=conv(xs,LP,'same');
    
    x3=medfilt1(xs,N);
    
    x4=wiener2(xs,[1,N]);
    
    %połączenia różnych filtrów dla lepszego wyniku
    x5=wiener2(medfilt1(wiener2(xs,[1,N]),N),[1,N]);
    x6=wiener2(wiener2(medfilt1(xs,N),[1,N]),[1,N]);
    x7=medfilt1(medfilt1(medfilt1(xs,N),N),N);

    blad(k,1)=ocena(x,x3);
    blad(k,2)=ocena(x,x4);
    blad(k,3)=ocena(x,x5);
    blad(k,4)=ocena(x,x7);
end

%plot(blad);
%legend;

%wychodzi na to że najlepszy jest:
%wiener2(medfilt1(wiener2(y,[1,N]),N),[1,N]) dla k=8

N=2*8+1;
W=wiener2(medfilt1(wiener2(xs,[1,N]),N),[1,N]);
L2=ocena(x,W);

FT=fftshift(fft(xs));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

subplot(211), plot(t,x,'b',t,W,'r');
subplot(212), plot(f,WA,'b');
L2