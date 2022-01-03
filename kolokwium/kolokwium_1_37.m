%% zad1
close all; clc; clear;

Fs=150;
t=0:1/Fs:6;


x=-1*((t>=0 & t<=1) | (t>=5 & t<=6)) + (1-4*(1-abs(t-3)/2)).*(abs(t-3)<2);

a0=-6/3;
FX=a0*ones(size(t))/2;

for n=1:76
   an=1/3*((3*(pi*n*sin((5*pi*n)/3)+6*cos((5*pi*n)/3)-12*cos(pi*n)-pi*n*sin((pi*n)/3)+6*cos((pi*n)/3)))/(pi^2*n^2)-(3*sin((pi*n)/3))/(pi*n)-(3*(sin(2*pi*n)-sin((5*pi*n)/3)))/(pi*n));
   bn=1/3*((3*(6*sin((5*pi*n)/3)-pi*n*cos((5*pi*n)/3)-12*sin(pi*n)+6*sin((pi*n)/3)+pi*n*cos((pi*n)/3)))/(pi^2*n^2)+(3*cos((pi*n)/3)-3)/(pi*n)+(3*(cos(2*pi*n)-cos((5*pi*n)/3)))/(pi*n));
   
   FX=FX+an*cos(n*pi*t/3)+bn*sin(n*pi*t/3);
   
   if n==10
       FX10=FX;
   end
   
end

plot(t,x,'.b',t,FX10,'k',t,FX,'g');

%% zad2
close all;clc;clear;
a=load('kolos/cor04.txt');
t=a(:,1)';
Fs=1/(t(2)-t(1));
x=a(:,2)';

tc=-57:1/Fs:57;

tp1=-3:1/Fs:3;
x1=0.75*exp(-(tp1).^2/(2*1*1));
xc1=xcorr(1-x.^3,1-x1.^3);
nr1=find(xc1==max(xc1(:)), 1, 'first');

tp2=-3:1/Fs:3;
x2=0.75*ones(size(tp2));
xc2=xcorr(x,x2);
nr2=find(xc2==max(xc2(:)), 1, 'first');

tp3=0:1/Fs:4;
x3=1*mod(tp3,4)/4;
xc3=xcorr(x.^5,x3.^5);
nr3=find(xc3==max(xc3(:)), 3, 'first');

tc(nr1)
tc(nr2)
tc(nr3)

subplot(211),plot(t,x,tp1,x1,'r',tp2,x2,'k',tp3,x3,'g');
subplot(212),plot(tc,xc3);

%% zad3
close all;clc;clear;
x=load('kolos/szum04.txt');
a=x(:,1)';
b=x(:,2)';
Fs=500;
t=linspace(0,length(a)/Fs,length(a));

ocena=@(x,y) (1/length(t)*sum(abs(x(:)-y(:))));
blad=zeros(49,4);

for N=3:2:101

%uśredniająca jednorodna (rzadko stosowana)
x1=conv(b,ones(1,N)/N,'same');

%uśredniająca ważona Gussem
N2=floor(N/2);
% 4 bo odchylenie nie powinno być większe niż 1/4 sygnału
LP=exp(-(-N2:N2).^2/(2*(N2/4).^2));
LP=LP/sum(LP(:));
x2=conv(b,LP,'same');

%mediana
x3=medfilt1(b,N);

%wiener
x4=wiener2(b,[1,N]);

blad(N2,1)=ocena(x1,a);
blad(N2,2)=ocena(x2,a);
blad(N2,3)=ocena(x3,a);
blad(N2,4)=ocena(x4,a);

%subplot(221), plot(t,x,'b',t,x1,'r');
%subplot(222), plot(t,x,'b',t,x2,'r');
%subplot(223), plot(t,x,'b',t,x3,'r');
%subplot(224), plot(t,x,'b',t,x4,'r');

end

FTa=fftshift(fft(a));
WAa=abs(FTa);

FT=fftshift(fft(b));
WA=abs(FT);
f=linspace(-Fs/2,Fs/2,length(t));

f0=5.8;
n=1;
Butterworth=(1.0)./(1+(f/f0).^(2*n));
c=real(ifft(ifftshift(Butterworth.*FT)));

L1=ocena(c,a);

subplot(211),plot(t,a,'b',t,c,'r');
subplot(212), plot(f,WAa,'r',f,WA,'b');
L1