clear all
clc;
Type=input('Enter type of filter 1:Lowpass 2:Highpass ');
display('for low pass Fp<Fs & for high pass Fp>Fs ');

Ap=input('Enter Passband attenuation ')
As=input('Enter Stopband attenuation ')
 Fp=input('Enter passband frequency ')
 Fs=input('Enter stopband frequency ')
 F=input('Enter sampling frequency ')

 
 fp=Fp/F
fst=Fs/F

switch Type
    case 1
        hs = fdesign.lowpass(fp,fst,Ap,As)
    case 2
        hs = fdesign.highpass(fst,fp,As,Ap)
    otherwise
        display('Error')
end

y2=butter(hs)

fvtool(y2)