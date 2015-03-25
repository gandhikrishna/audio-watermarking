m = 1/sqrt(2.*pi);   % max. value for normal function
[y, Fs] = wavread('Cena.wav'); % read .wav file
x = 0;
len = length(y);
ns = floor(len/200);

for j=0:199
    
    temp1 = (ns*j)+1;
    temp2 = ns*(j+1);
    
    smpl = y(temp1:temp2);    %divide signal into 200 frames

    dct_smpl = dct(smpl);     %take dft of each frame
%    mag = abs(dft_smpl);      %magnitude response
 %   mag_n = mag;
  %  ph = angle(dft_smpl);     %Phase response

    [c,k] = max(dct_smpl);         %detection of peak value
%    if k == 1
 %       mag_n(1) = 0;
  %      [c1,k1] = max(mag_n);
   %     mag_n(1) = c;
    %    c = c1;
     %   k = k1;
   % end
    
    x(j+1) = m*rand;          %generate watermark stream
    v = dct_smpl(k);
    alpha = 0.1;
    vn = v*(1+alpha*x(j+1));  %new weighted value

    dct_smpl(k) = vn;
 %   mag_n(ns+2-k) = vn;

  %  fn1 = mag_n.*(cos(ph)+(i.*sin(ph)));  %generate new dft

    smpl_m = (idct(dct_smpl));                 %get watermarked frame by inverse fft
    smpl_n(temp1:temp2) = (smpl_m);

end
subplot(2,1,1)
plot(y)
subplot(2,1,2)
plot(smpl_n)

%smpl_final = awgn(smpl_final,35);
wavwrite(smpl_n,Fs,'wm_try2.wav');

sm = 0;
sm1 = 0;
for vr = 1:length(smpl_n)
    sm = sm + y(vr)^2;
    sm1 = sm1 + (y(vr) - smpl_n(vr))^2;
end

snr = 10*log10(sm/sm1);
display('Signal to noise ratio (SNR) is :');
display(snr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[z, Fk] = wavread('wm_try2.wav');
[y, Fs] = wavread('Cena.wav');

len = length(y);
ns = floor(len/200);

for j=0:199

    temp1 = (ns*j)+1;
    temp2 = ns*(j+1);
    
    smpl_z = z(temp1:temp2);
    smpl = y(temp1:temp2);


    dct_z = dct(smpl_z);
%    mag_z = abs(dft_z);
%    fn_t = mag_z;
 %   ph_z = angle(dft_z);

    f1 = dct(smpl);
%    f = abs(f1);
%    fn = f;
%    p1 = angle(f1);

    [c,k] = max(f1);
 %   if k == 1
 %       f(1) = 0;
 %       [c2,k2] = max(f);
 %       f(1) = c;
 %       c = c2;
 %       k = k2;
  %  end
    
    [c_z,k_z] = max(dct_z);
 %   if k_z == 1
  %      mag_z(1) = 0;
   %     [c3,k3] = max(mag_z);
    %    mag_z(1) = c;         
     %   c_z = c3;
      %  k_z = k3;
 %   end

    z_t(j+1) = c-c_z;

    alpha = 0.1;
    x_z(j+1) = ((c_z/c)-1)/alpha;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(x);
p1=N*sum(x.*x_z)-sum(x)*sum(x_z);
p2=N*sum(x.*x) - (sum(x))^2;
p3=N*sum(x_z.*x_z) - (sum(x_z))^2;
p4=sqrt(p2*p3);
corelation=p1/p4;
display('Correlation is :'); 
display(corelation);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%crp1 = input('Input Crop percentage :');
%crp = 100/crp1;
%lx = floor(len/crp);
y1=downsample(z,2);
z=upsample(y1,2);
%smpl_final(1:lx) = y(1:lx);
smpl_new = z;
wavwrite(smpl_new,Fs,'wm_try3.wav');
ns = floor(len/200);

for j=0:199

    temp1 = (ns*j)+1;
    temp2 = ns*(j+1);
    
    smpl_z = smpl_new(temp1:temp2);
    smpl = y(temp1:temp2);


    dct_z_new = dct(smpl_z);
  %  mag_z = abs(dft_z);
  %  fn_t = mag_z;
  %  ph_z = angle(dft_z);

    f1 = dct(smpl);
 %   f = abs(f1);
 %   fn = f;
 %   p1 = angle(f1);

    [c,k] = max(f1);
    [c_z,k_z] = max(dct_z_new);

    z_t(j+1) = c-c_z;

    alpha = 0.1;
    x_z(j+1) = ((c_z/c)-1)/alpha;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% New Correlation

N=length(x);
p1=N*sum(x.*x_z)-sum(x)*sum(x_z);
p2=N*sum(x.*x) - (sum(x))^2;
p3=N*sum(x_z.*x_z) - (sum(x_z))^2;
p4=sqrt(p2*p3);
corelation=p1/p4;
display('The new correlation after atack is :')
display(corelation);