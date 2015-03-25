clear all
close all

Ap=input('Enter passband attenuation=');
As=input('Enter stopband attenuation=');
Fp=input('Enter passband frequency=');
Fs=input('Enter stopband frequency=');
F=input('Enter sampling frequency=');
T=1/F

wp=2*pi*Fp
ws=2*pi*Fs

ohmp=wp/T
ohms=ws/T

x1=10^(0.1*As)-1
x2=10^(0.1*Ap)-1
x3=log(x1/x2)
x4=log(ohmp/ohms)
N=ceil(x3/(2*x4))

x5=x2^(1/(2*N))
ohmc=ohmp/x5

[num,den]=maxflat(0,N,ohms)
%[num1,den1]=lp2hp(num,den,ohms)
fvtool(den,num)