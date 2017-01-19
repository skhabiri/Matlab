clear all
close all
format long

% Modulator Parameters
alpha1 = 1/16; %%Critical!
alpha2 = 1/64; %%Critical!
beta1  = 1/32;
beta2  = 1/16;
beta3  = 1/8;
beta4  = 1/4;
beta5  = 1/2;

len       = 2^16;      % Total number of points
fs        = 48000;     % Sampling frequency
overs     = 64;       % Oversampling ratio
freq      = 1*1031.25;   % 5859.375;  %1031.25;   % 22*48e3*64/2^16=1031.25
Outpower=2.5; % maximum power is PVdd^2/16
ClassdVdd=3.6;
PVdd=8; %boosted power stage
amplitude = (sqrt(Outpower*16)/(3*2.5))      % Sinewave Amplitude


% Modeling the input sinewave

  t  = (0:1:len-1)/(fs*overs);          % Time scale
  x  = amplitude*sin(2*pi*freq*t);      % Input sinewave
  xz = fft(x);	             	        % Input FFT

% Initializing output of modulator , output vector and internal states of the modulator

  y = zeros(1,len);                    % output of the modulator
  i = zeros(1,5);		       % new values
  i_o = zeros(1,5);                    % old values
  out = 0;
  
% Class-D Params
UGB_loop=18e4;%15e4;

m=2; % number of DAC output bit;


C=10e-12;
Rf2ndArray=[0.8*1/(2*pi*(UGB_loop)*C) 1*1/(2*pi*(UGB_loop)*C) 1.2*1/(2*pi*(UGB_loop)*C)];
Rf2nd=Rf2ndArray(2);

Rf=1.5*Rf2nd/0.875;
%Rf=180e3;
C1p=2*C;

Rf2=Rf2nd*2*1.8;
%Rf2=340e3;
C2p=C/2;

Rf3=Rf2nd*4;
%Rf3=420e3;
C3=C/4;

Rin2nd=5*Rf/3;
Rin0=3*Rin2nd
Rin1=3*Rin2nd/2

fac=1e3;

%fsPWM=0.125*(fs*overs);
fsPWM=0.125*(fs*overs);
RdsonP=0.15;
RdsonN=0.15;
Rmetal=0.3;
wu1=17e6*2*pi;
p3db1=260*2*pi;
g01=wu1/p3db1;

wu2=10e6*2*pi;
p3db2=260*2*pi;
g02=wu2/p3db2;

wu3=10e6*2*pi;
p3db3=520*2*pi;
g03=wu3/p3db3;

Ts=0.05/(fs*overs);


% Modeling the Fifth Order Digital Sigma-Delta Modulator 
%{
   for n=2:1:len, 

      i(1) = (beta1*x(n-1) + i(1) - beta1*out);
      i(2) = (beta2*i_o(1) + i(2) - beta2*out - beta2*i(3)*alpha1);
      i(3) = (beta3*i_o(2) + i(3) - beta3*out);
      i(4) = (beta4*i_o(3) + i(4) - beta4*out - beta4*i(5)*alpha2);
      i(5) = (beta5*i_o(4) + i(5) - beta5*out);
      i5(n)=i(5);

      % Modeling the output comparator
      if( i(5) >= 0 )
        out =  1;
      else
        out = -1;
      end
 
      y(n) = out;
 
      % Updating old values	
      i_o(1) = i(1);
      i_o(2) = i(2);
      i_o(3) = i(3);
      i_o(4) = i(4);
      i_o(5) = i(5);
 
   end              
%}
%    ----------------------------- end of the loop ---------------------
 sim('Sim_mbitDACPlusClassD')
 y=y';
   [max_samp_fft,i] = max(abs(xz(1:len/4)));      % Calculating where the sine is
   snr = SNR1(y,i);                               % Calculating SNR, result in snr vector
 %  xxx = (fs/2)*len/(overs*fs);                     % xxx is the bin number for 24 KHz (fs/2)
   xxx = round((20000)*len/(overs*fs));                     % xxx is the bin number for 24 KHz (fs/2)
   
   % Plotting results
   gg = [7500 75000]; %5th order slope reference line (100dB/dec)
   hh = [-150 -50];
   
   
   figure
   plot((2:xxx)*fs*overs/len,snr(2:xxx))
   grid
   ylabel('Output Signal-to-Noise Ratio <dB>')
   xlabel('Frequency <Hz>')
   title('Output Signal-to-Noise Ratio vs. Frequency for the 5th Order DSDM')
   axis([0 20000 -160 -80]);

   

  
 
  
   figure
   f = 1:len;
   f = f/len*fs*overs;
   argu = abs(fft(y.*blackman(len)'));
   max_argu = max(argu);
   mm = 20*log10(argu/max_argu);
   semilogx(f(1:len/2),mm(1:len/2))
   hold 
   semilogx(gg,hh,'r')
   grid
   xlabel('Frequency <Hz>')
   ylabel('Output FFT <dB>')
   text(3000,-80,'100 dB/dec ->')
   title('5th Order DSDM Output FFT')

   figure
   plot(f(1:xxx),mm(1:xxx),f(1:xxx),snr(1:xxx))
   axis([0 20000 -150 6]);
   grid
   xlabel('Frequency <Hz>')
   ylabel('Output FFT <dB>')
   title('Zoom-In 5th Order DSDM Output FFT')

   % Printing the SNR Value at 24 KHz
  SNR_in=snr(xxx)
  
   M=2^20;
   win=kaiser(M,20);
   %Ts=0.05/(fs*overs);
   Fs=1/Ts;
   fin=freq;
  
   outwin=simout1(100000:M+100000-1).*win;
   outfft=fft(outwin);
   outabs=(outfft.*conj(outfft));
   outdb=10*log10(outabs/max(outabs));
   
   figure
   f1=0:Fs/M:Fs-Fs/M;
   plot(f1(1:xxx),outdb(1:xxx))
   axis([0 20000 -150 6]);
   grid
   xlabel('Frequency <Hz>')
   ylabel('Class-D Output FFT <dB>')
   title('Zoom-In 5th Order DSDM plus Class-D Output FFT')
  
   S=sum(outabs(max(8,ceil(fin/Fs*M))-7:ceil(fin/Fs*M)+7));
   DPN=sum(outabs(1:ceil(2e4/Fs*M)))-S;
   THD_ClassD_Out=sqrt(DPN/S)*100
   
   snr_classD = SNR1(y_classD',i);
   SNDR_out=snr_classD(xxx)
   
   snr_classD_Aw = SNR1(y_classD_Aw',i);
   SNDR_out_Aw=snr_classD_Aw(xxx)
   
   SNR_at_3W=SNDR_out - 10*log10(3/(Outpower))
   SNR_at_3W_Aw=SNDR_out_Aw - 10*log10(3/(Outpower))
  
  
