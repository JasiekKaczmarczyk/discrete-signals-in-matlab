%%
close all;clc;clear;

Fs=100;
t=0:1/Fs:6;

x=1*((t>=0 & t<=1) | (t>=5 & t<=6)) + (4*(1-abs(t-3)/2)-2).*(abs(t-3)<2);

a0=2/3;
FX=a0*ones(size(t))/2;

for n=1:70
   an=1/3*(-(6*(pi*n*sin((5*pi*n)/3)+3*cos((5*pi*n)/3)-6*cos(pi*n)-pi*n*sin((pi*n)/3)+3*cos((pi*n)/3)))/(pi^2*n^2)+(3*sin((pi*n)/3))/(pi*n)+(3*(sin(2*pi*n)-sin((5*pi*n)/3)))/(pi*n));
   bn=1/3*(-(6*(3*sin((5*pi*n)/3)-pi*n*cos((5*pi*n)/3)-6*sin(pi*n)+3*sin((pi*n)/3)+pi*n*cos((pi*n)/3)))/(pi^2*n^2)-(3*cos((pi*n)/3)-3)/(pi*n)-(3*(cos(2*pi*n)-cos((5*pi*n)/3)))/(pi*n));
   
   FX=FX+an*cos(n*pi*t/3)+bn*sin(n*pi*t/3);
   
   if n==9
      FX9=FX; 
   end
end

plot(t,x,'.b',t,FX9,'r',t,FX,'g');

%%
close all;clc;clear;
a=load('cor08.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-110:1/Fs:110;

tp1=-3:1/Fs:3;
x1=(1-abs(tp1)/3);
xc1=xcorr(x,x1)+xcorr(1-x,1-x1);
nr1=find(xc1==max(xc1(:)), 1, 'first');

tc(nr1)

tp2=0:1/Fs:5;
x2=mod(tp2,5)/5;
xc2=xcorr(x.^3,x2.^3);
nr2=find(xc2==max(xc2(:)), 3, 'first');

%wypisanie czasu największej korelacji
tc(nr2)


tp3=-5:1/Fs:5;
x3=0.65*ones(size(tp3));
xc3=xcorr(x,x3);
nr3=find(xc3>0.9*max(xc3(:)), 1, 'first');

tc(nr3)

subplot(211),plot(t,x,tp1,x1,'g',tp2,x2,'r',tp3,x3,'y');
subplot(212),plot(tc,xc3);

%% zad3
close all;clc;clear;
x=load('szum07.txt');
a=x(:,1)';
b=x(:,2)';
Fs=200;
t=linspace(0,length(a)/Fs,length(a));

ocena=@(x,y) (1/length(t)*sum(abs(x(:)-y(:))));

FTa=fftshift(fft(a));
WAa=abs(FTa);

FT=fftshift(fft(b));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

f0=3.9;
n=1;
Butterworth=(1.0)./(1+(f/f0).^(2*n)).*(abs(f)<=10);
c=real(ifft(ifftshift(Butterworth.*FT)));


L1=ocena(a,c);

subplot(211),plot(t,a,'b',t,c,'r');
subplot(212), plot(f,WAa,'r',f,WA,'b',f,100000*Butterworth,'g');
L1

%% zad3
close all;clc;clear;
x=load('szum07.txt');
a=x(:,1)';
b=x(:,2)';
Fs=200;
t=linspace(0,length(a)/Fs,length(a));

ocena=@(x,y) (1/length(t)*sum(abs(x(:)-y(:))));

FTa=fftshift(fft(a));
WAa=abs(FTa);

FT=fftshift(fft(b));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

BStop=1.0*(abs(f)<1);
FT_n=BStop.*FT;
c=real(ifft(ifftshift(FT_n)));
L1=ocena(a,c);

subplot(211),plot(t,a,'b',t,c,'r');
subplot(212), plot(f,WAa,'r',f,WA,'b',f,100000*BStop,'g');
L1