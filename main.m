
 clc; clear; close all;

  
 M = 16; % order of the modulation
 n = 16; % number of bits in bitstream
 k=log2(M); % bits per symbol
 EbNo = 10; % to calculate snr
 fc=100; %frequency of the carrier
 Tb=1/n;		        % Bit duration time [s]
 

fs=1000;     	    % sampling frequency [Hz]
fn=fs/2;            % Nyquist frequency [Hz]
Ts=1/fs;	        % Sampling time [s]
t=[0:Ts:(n*Tb)-Ts]';% Time vector initialization
 
 
 
% fsym=1/Tb; %freq of bits;% rolloff = 0.5;
% rolloff = 0.5;
% snr = 20;
% kanal = 2;
% BR = 64000;


stream=RandBitStream(n); %generating random bitstream


% Modulator

%mapping
[ map, complex_constell ]=gray(M); %constellation


grouped_bits=bitgrp(stream, M);
symbols=binary2dec(grouped_bits);
mapped=map2gray(symbols, map);

%splitting
[I, Q]=split_stream(mapped);


% ??
fig = scatterplot(complex_constell);
hold on;
pause;
scatterplot(mapped(:,2:3), 1, 0,'ro',fig);
pause;




% creating modulating signals sI and sQ

len=n/k; %length of I or Q in symbols
rep=fs/len; %no of repetitions of a single bit


x=1:1:(len+1)*(1/Ts*k);
for i=1:len
    for j=1:.1:i+1;
        sI(x(i*rep:(i+1)*rep))=I(i);
        sQ(x(i*rep:(i+1)*rep))=Q(i);
    end
end
sI=sI(rep:end-1);
sQ=sQ(rep:end-1);




hold off;
plot(t, sI, 'r'); ylabel ('I amplitude');  xlabel ('t[s]');
pause
plot(t, sQ, 'g'); ylabel ('Q amplitude'); xlabel ('t[s]');
pause

% Carriers
carrier_I=cos(2*pi*fc*t);
carrier_Q=-sin(2*pi*fc*t);

% Modulating
modulated_I=sI'.*carrier_I;
modulated_Q=sQ'.*carrier_Q;

% Summing I and Q
output_signal=modulated_I+modulated_Q;

plot(t, modulated_I); title ('Modulated I'); ylabel ('I amplitude');  xlabel ('t[s]');
pause;
plot(t, modulated_Q); title ('Modulated Q'); ylabel ('Q amplitude');  xlabel ('t[s]');
pause
plot(t, output_signal); title ('Output Signal - sum of modulated I and Q'); ylabel ('Amplitude');  xlabel ('t[s]');
pause;


% %transmittiing through AWGN

snr = EbNo + 10*log10(k);


% receivedSignal = awgn(mapped(:,2:3), snr, 'measured'); %% to be replaced
%
% output_signal + noise !!!!!!!!!!!!!!!!
%noise=randn(size(mapped(:,2:3))); % random noise generation
%constant=std(mapped(:,2:3))/(std(noise)*10^(snr/20));
%receivedSignal=mapped(:,2:3) + noise*constant; %output of transmitter
%noise1=noise*constant;


%demodulation

I_recovered=output_signal.*carrier_I;
Q_recovered=output_signal.*carrier_Q;

%before LPF
plot(t, I_recovered); title ('Recovered I'); ylabel ('I amplitude');  xlabel ('t[s]');
pause


plot_fft(I_recovered, fs, 'r');


%filtration

filtered_I=lpf(I_recovered, fc, fs);
filtered_Q=lpf(Q_recovered, fc, fs);


filtered_I=2*filtered_I;
filtered_Q=2*filtered_Q;

%averaging

sym_num=n/k;
demodulated_I=mean(reshape(filtered_I, fs/sym_num, []))';
demodulated_Q=mean(reshape(filtered_Q, fs/sym_num, []))';
receivedSignal=horzcat(demodulated_I, demodulated_Q)
mappedSignal=mapped(:,2:3)

%scatterplot(receivedSignal, 1, 0, 'g.', fig);

recovered_constellation=rec_constell(receivedSignal, map);

%scatterplot((2.*round((recovered_constellation+1)/2)-1), 1, 0, 'ko', fig);

recovered_symbols=gray2symbols(recovered_constellation, map);
recovered_bits=symbols2bits(recovered_symbols);

[bit_errors, ber]=ber_calc(recovered_bits, stream);



%theoretical ber, some approximations taken // wiki
x=sqrt(3*k*EbNo/(M-1));
theoretical_ber=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));


bit_errors
ber
theoretical_ber
length(recovered_bits)
