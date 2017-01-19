close all
clear all
sample_rate = 48000;
len = 1024;	%number of samples
freqrange = (1:len)*sample_rate/len;
a=round(len*20000/sample_rate);	%bin# representing 20KHz

load -ascii left.txt
fs = 48000;
start_offset = 100;
number_of_bits = 20;
left_output = left(start_offset:len+start_offset);

for n=1:len
 if (left_output(n) >= 2^(number_of_bits-1))
  left_output(n) = left_output(n)-2^number_of_bits ; %decimal negative numbers
 end
end

data = left_output;
argu = 20*log10(abs(fft(data(1:len).*blackman(len))));
max_argu = max(argu);
[max_samp_fft,ii] = max(argu(1:a)); % Calculating where the fundamental sine is


%snr = SNR3_aw(data(1:len)',ii,sample_rate);


[temp,fftsize] = size(data(1:len)');
no = zeros(1,fftsize);
fso2 = fftsize/2-1;
fundwid = 10;	% signal window roughly
fundwidh = 10;	% harmonics window roughly

sigmin = ii-fundwid;
if(sigmin < 0)
	sigmin = 1;
end

sigmax = ii+fundwid;

if(sigmax>fso2)
	sigmax = fso2;
end

noisefloor = 1e-40;
freq = ((1:fftsize)/fftsize).*sample_rate;

% This is in dB:
awlog = (-(log(freq/2560).^2.25)*2.1) + 1.3;

% So convert to linear:
aw = abs(10.^(awlog./20));

% And multiply the linear data
xy = fft(data(1:len)'.*blackman(fftsize)');

xy = xy.*aw;
%% End of A-weighting section


noiso = noisefloor;

outdata = (abs(xy).^2);

for i =sigmin:sigmax
	sigpwr = noiso + outdata(i);	%This is signal power
end

no(1)=noisefloor;
signoise = (outdata((ii)-fundwid)+outdata((ii)+fundwid))/2;
if(((3*ii)+fundwidh) < fso2)
no3rd = (outdata((3*ii)-fundwidh) + outdata((3*ii)+fundwidh))/2;
end
if(((2*ii)+fundwidh) < fso2)
no2nd = (outdata((2*ii)-fundwidh) + outdata((2*ii)+fundwidh))/2;
end
if(((5*ii)+fundwidh) < fso2)
no5th = (outdata((5*ii)-fundwidh) + outdata((5*ii)+fundwidh))/2;
end
if(((7*ii)+fundwidh) < fso2)
no7th = (outdata((7*ii)-fundwidh) + outdata((7*ii)+fundwidh))/2;
end
if(((9*ii)+fundwidh) < fso2)
no9th = (outdata((9*ii)-fundwidh) + outdata((9*ii)+fundwidh))/2;
end

for i=10:fso2	%skip DC bins; integration
	if(i<(ii+fundwid) & i>(ii-fundwid))
		no(i) = no(i-1)+signoise;
        elseif (i<((3*ii)+fundwidh) & i>((3*ii)-fundwidh))
                no(i) = no(i-1) + no3rd;
        elseif (i<((2*ii)+fundwidh) & i>((2*ii)-fundwidh))
                no(i) = no(i-1) + no2nd;
%        elseif (i<((5*ii)+fundwidh) & i>((5*ii)-fundwidh))
%                 no(i) = no(i-1) + no5th;
%         elseif (i<((7*ii)+fundwidh) & i>((7*ii)-fundwidh))
%                 no(i) = no(i-1) + no7th;
%         elseif (i<((9*ii)+fundwidh) & i>((9*ii)-fundwidh))
%                 no(i) = no(i-1) + no9th;
	else
		no(i) = no(i-1) + outdata(i);
	end  
end

snr = 10*log10(no/sigpwr); %it's not taking the harmonic spikes into calculation, only edges of harmonics are considered



figure
plot(freqrange(1:a),argu(1:a)-max_samp_fft,freqrange(1:a),snr(1:a));
left_snr = -snr(a)+7.7
axis([0 20000 -150 0]);
grid
xlabel('Frequency <Hz>')
ylabel('Output FFT <dB>')
title(['20KHz single ended ADC Decimation Filter Output',', SNR=',num2str(left_snr),' (FSdB)'])

snr_file_id = fopen('adc_snr_results.mat','w');
fprintf(snr_file_id,'%d\n',-1*round(left_snr));
fclose(snr_file_id)
