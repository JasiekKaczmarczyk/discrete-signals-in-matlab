%% zad 1
close all; clc; clear;
Fs=150;
t=0:1/Fs:6;


x=(-1*(t>=0 & t<=2)) + (-1*(t>=4 & t<=6)) + (3*(1-abs(t-3))-2).*(abs(t-3)<1);

a0=2/3;
FS=a0*ones(size(t))/2;

for n=1:120 
  an=1/3*(-(3*(2*pi*n*(4*sin((4*pi*n)/3)-sin((2*pi*n)/3))+9*cos((4*pi*n)/3)-9*cos((2*pi*n)/3)))/(pi^2*n^2)-(3*sin((2*pi*n)/3))/(pi*n)-(3*(sin(2*pi*n)-sin((4*pi*n)/3)))/(pi*n));
  bn=1/3*(-(27*sin((4*pi*n)/3)-24*pi*n*cos((4*pi*n)/3)-27*sin((2*pi*n)/3)+6*pi*n*cos((2*pi*n)/3))/(pi^2*n^2)+(3*cos((2*pi*n)/3)-3)/(pi*n)+(3*(cos(2*pi*n)-cos((4*pi*n)/3)))/(pi*n));
  
  FS=FS+an*cos(n*pi*t/3)+bn*sin(n*pi*t/3);
  
  if n==8
      x8=FS;
  end
 
end

plot(t,x,'y.',t,x8,'b',t,FS,'black');

%% zad cn
close all; clc; clear;
Fs=150;
t=0:1/Fs:6;


x=(-1*(t>=0 & t<=2)) + (-1*(t>=4 & t<=6)) + (3*(1-abs(t-3))-2).*(abs(t-3)<1);

c0=2/3;
FS=c0;

for n=1:120
  cn=1/12*(-12*((pi*n+9j)*sin(2*pi*n/3)+(9-1j*pi*n)*cos(2*pi*n/3)-18j*sin(pi*n/2)-18*cos(pi*n/2)+(9j-pi*n)*sin(pi*n/3)+(1j*pi*n+9)*cos(pi*n/3))/(pi*pi*n*n) -6*(sin(pi*n)-1j*cos(pi*n)-sin(2*pi*n/3)+1j*cos(2*pi*n/3)+sin(pi*n/3)-1j*cos(pi*n/3)+1j)/(pi*n));
  
  FS=FS+cn*exp(1j*n*pi*t/6);
  
  if n==8
      x8=FS;
  end
  
  if n==120
     x120=FS; 
  end
end

plot(t,x,'y.',t,x8,'b',t,x120,'black');

%%
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

%subplot(211), plot(t,x,tp3,x3);
%subplot(212), plot(tc,xc3);

%%
close all; clc;clear;

a=load('kolos/szum08.txt');
x0=a(:,1)';
x=a(:,2)';

Fs=400;
t=linspace(0,1/Fs*length(x),length(x));

ocena=@(x,y) ((1/length(t))*sum(abs(x(:)-y(:))));
blad=zeros(25,4);

for N=3:2:51

%uśredniająca jednorodna (rzadko stosowana)
x1=conv(x,ones(1,N)/N,'same');

%uśredniająca ważona Gussem
N2=floor(N/2);
% 4 bo odchylenie nie powinno być większe niż 1/4 sygnału
LP=exp(-(-N2:N2).^2/(2*(N2/4).^2));
LP=LP/sum(LP(:));
x2=conv(x,LP,'same');

%mediana
x3=medfilt1(x,N);

%wiener
x4=wiener2(x,[1,N]);

blad(N2,1)=ocena(x0,x1);
blad(N2,2)=ocena(x0,x2);
blad(N2,3)=ocena(x0,x3);
blad(N2,4)=ocena(x0,x4);

%subplot(221), plot(t,x0,'b',t,x1,'r');
%subplot(222), plot(t,x0,'b',t,x2,'r');
%subplot(223), plot(t,x0,'b',t,x3,'r');
%subplot(224), plot(t,x0,'b',t,x4,'r');

end

N=53;
%uśredniająca ważona Gussem
N2=floor(N/2);
LP=exp(-(-N2:N2).^2/(2*(N2/4).^2));
LP=LP/sum(LP(:));
x_naj=conv(x,LP,'same');
ocena(x0,x_naj)

FT=fftshift(fft(x0));
WA=abs(FT);
f=linspace(-Fs,Fs,length(t));

plot(blad);

%subplot(211), plot(t,x,'b',t,x_naj,'r');
%subplot(212), plot(f,WA,'b');