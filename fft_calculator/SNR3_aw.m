function   snr = SNR3_aw(y,frq,fs)

[temp,fftsize] = size(y);
no = zeros(1,fftsize);
nod = zeros(1,fftsize);

fso2 = fftsize/2-1;
fundwid = 10;	% signal window roughly
fundwidh = 10;	% harmonics window roughly

sigmin = frq-fundwid;
if(sigmin < 0)
	sigmin = 1;
end

sigmax = frq+fundwid;

if(sigmax>fso2)
	sigmax = fso2;
end

noisefloor = 1e-40;
freq = ((1:fftsize)/fftsize).*fs;

% This is in dB:
awlog = (-(log(freq/2560).^2.25)*2.1) + 1.3;

% So convert to linear:
aw = abs(10.^(awlog./20));

% And multiply the linear data
xy = fft(y.*blackman(fftsize)');

xy = xy.*aw;
%% End of A-weighting section


sigpwr = noisefloor;

outdata = (abs(xy).^2);

for i =sigmin:sigmax
	sigpwr = sigpwr + outdata(i);	%This is signal power
end
no2nd=0;
no3rd=0;
no5th=0;
no7th=0;
no9th=0;

nod(1)=noisefloor;
no(1)=noisefloor;

signoise = (outdata((frq)-fundwid)+outdata((frq)+fundwid))/2;
if(((3*frq)+fundwidh) < fso2)
no3rd = (outdata((3*frq)-fundwidh) + outdata((3*frq)+fundwidh))/2;
end
if(((2*frq)+fundwidh) < fso2)
no2nd = (outdata((2*frq)-fundwidh) + outdata((2*frq)+fundwidh))/2;
end
if(((5*frq)+fundwidh) < fso2)
no5th = (outdata((5*frq)-fundwidh) + outdata((5*frq)+fundwidh))/2;
end
if(((7*frq)+fundwidh) < fso2)
no7th = (outdata((7*frq)-fundwidh) + outdata((7*frq)+fundwidh))/2;
end
if(((9*frq)+fundwidh) < fso2)
no9th = (outdata((9*frq)-fundwidh) + outdata((9*frq)+fundwidh))/2;
end

for i=10:fso2	%skip DC bins; integration
	if(i<(frq+fundwid) & i>(frq-fundwid))
		no(i) = no(i-1)+signoise;
        elseif (i<((3*frq)+fundwidh) & i>((3*frq)-fundwidh))
                no(i) = no(i-1) + no3rd;
        elseif (i<((2*frq)+fundwidh) & i>((2*frq)-fundwidh))
                no(i) = no(i-1) + no2nd;
        elseif (i<((5*frq)+fundwidh) & i>((5*frq)-fundwidh))
                 no(i) = no(i-1) + no5th;
         elseif (i<((7*frq)+fundwidh) & i>((7*frq)-fundwidh))
                 no(i) = no(i-1) + no7th;
         elseif (i<((9*frq)+fundwidh) & i>((9*frq)-fundwidh))
                 no(i) = no(i-1) + no9th;
	else
		no(i) = no(i-1) + outdata(i);
	end  
end

snr = 10*log10(no/sigpwr); %it's not taking the harmonic spikes into calculation, only edges of harmonics are considered, sigpwr is scalar, but no is vector

for i=10:fso2	%skip DC bins; integration
	if(i<(frq+fundwid) & i>(frq-fundwid))
		nod(i) = nod(i-1)+signoise;
	else
		nod(i) = nod(i-1) + outdata(i);
	end  
end
sndr = 10*log10(nod/sigpwr); %it's not taking the harmonic spikes into calculation, only edges of harmonics are considered
