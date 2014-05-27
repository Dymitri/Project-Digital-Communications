function [jml_bit_err x]=qam(M,rolloff,snr,kanal,BR)
% clc; clear; close all;
% 
% M = 256;
% rolloff = 0.5;
% snr = 20;
% kanal = 2;
% BR = 64000;

nsamp = 4; % Oversampling rate
Tb = 1/BR;     %perioda 1 bit
tau = 1.0e-6;
shift = round(tau/Tb);

% Modulator
% M = 16;
x = randint(1,BR*0.02,M);
y = qammod(x,M);

% Filter Definition
% Define filter-related parameters.
filtorder = 40; % Filter order
delay = filtorder/(nsamp*2); % Group delay (# of input samples)
% rolloff = 0.25; % Rolloff factor of filter

% Create a square root raised cosine filter.
rrcfilter = rcosine(1,nsamp,'fir/sqrt',rolloff,delay);

% Transmitted Signal
% Upsample and apply square root raised cosine filter.
ytx = rcosflt(y,1,nsamp,'filter',rrcfilter);


% ==== Kanal ==============================================================
if kanal==1
    ytx = awgn(ytx,snr,'measured');                   
elseif kanal==2
    ytx=ytx';
    c1 = rayleighchan;                    %Kanal Rayleigh tap 1
    reset(c1,[12345;54321]);
    c1.ResetBeforeFiltering = 0;
    c2 = rayleighchan;                    %Kanal Rayleigh tap 2
    reset(c2,[54321;12345]);
    c2.ResetBeforeFiltering = 0;
                
    fad1 = filter(c1,ones(1,length(ytx))); fad1 = abs(fad1);
    fad2 = filter(c2,ones(1,length(ytx))); fad2 = abs(fad2);
    faded_sig1 = fad1.*ytx;
    faded_sig2 = fad2.*[zeros(1,shift) ytx(1:end-shift)];

    sigPower1 = sum(abs(faded_sig1(:)).^2)/length(faded_sig1(:));
    sigPower1db = 10*log10(sigPower1);
    sigPower2 = sum(abs(faded_sig2(:)).^2)/length(faded_sig2(:));
    sigPower2db = 10*log10(sigPower2);
    att = (10^(-3/10))*sigPower1/sigPower2;
    faded_sig2 = faded_sig2.*sqrt(att);

    faded_Sig = faded_sig1 + faded_sig2;

    ytx=awgn(faded_Sig,snr,'measured');
    ytx=ytx';
end

% Received Signal
% Filter received signal using square root raised cosine filter.
yrx = rcosflt(ytx,1,nsamp,'Fs/filter',rrcfilter);
yrx = downsample(yrx,nsamp); % Downsample.
yrx = yrx(2*delay+1:end-2*delay); % Account for delay.

% Demodulate to recover the message.
z = qamdemod(yrx,M);

% Check symbol error rate.
[num,rt]= symerr(x,z');

jml_bit_err = num;
ser = rt;