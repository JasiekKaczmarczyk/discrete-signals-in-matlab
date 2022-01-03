%%
close all; clc; clear;

Fs=125;
t=0:1/Fs:4;

x=-1*(t>=0 & t<=1) + (-4*t+3).*(t>=1 & t<=2) + (4*t-11).*(t>=2 & t<=3) -1*(t>=3 & t<=4);

plot(t,x,'r')