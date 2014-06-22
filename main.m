%               QAM demo
%Authors:
%    Przemys³aw Dymitrowski
%    Piotr Szkotak
%    Electronics and Telecommunications III year
%    06.2014
%
%   To use it adjust parameters below or use dedicated Graphical User
%   Interface by executing gui.m file
%
%
%


clc; clear; close all;

  
 M = 16; % order of the modulation
 n = 16; % number of bits in bitstream
 k=log2(M); % bits per symbol
 EbNo = 20; % to calculate snr
 fc=40; %frequency of the carrier
 Tb=1/n;		        % Bit duration time [s]
 

fs=500;     	    % sampling frequency [Hz]
fn=fs/2;            % Nyquist frequency [Hz]
Ts=1/fs;	        % Sampling time [s]
t=[0:Ts:(n*Tb)-Ts]';% Time vector initialization
 
stream=RandBitStream(n); %generating random bitstream


% Modulator

%mapping
[ map, complex_constell ]=gray(M); %creating constellation
grouped_bits=bitgrp(stream, M); %grouping bits
symbols=binary2dec(grouped_bits); % creating symbols
mapped=map2gray(symbols, map); %mapping symbols to the constellation points

%splitting
[I, Q]=split_stream(mapped);



% creating modulating signals sI and sQ

len=n/k; %length of I or Q in symbols
rep=fs/len; %number of repetitions of a single bit (to obtain square wave)

x=1:1:(len+1)*(1/Ts*k);
for i=1:len
    for j=1:.1:i+1;
        sI(x(i*rep:(i+1)*rep))=I(i);
        sQ(x(i*rep:(i+1)*rep))=Q(i);
    end
end
sI=sI(rep:end-1); %modulating signals
sQ=sQ(rep:end-1);



% Carriers
carrier_I=cos(2*pi*fc*t);
carrier_Q=-sin(2*pi*fc*t);

% Modulating
modulated_I=sI'.*carrier_I;
modulated_Q=sQ'.*carrier_Q;

% Summing I and Q
output_signal=modulated_I+modulated_Q;

%plots 
plot(t, sI, 'r'); title ('Modulating signal I'); ylabel ('I amplitude');  xlabel ('t[s]'); pause;
plot(t, modulated_I, 'r'); title ('Modulated I'); ylabel ('I amplitude');  xlabel ('t[s]'); pause;

plot(t, sQ, 'g'); title ('Modulating signal Q'); ylabel ('Q amplitude'); xlabel ('t[s]'); pause;
plot(t, modulated_Q, 'g'); title ('Modulated Q'); ylabel ('Q amplitude');  xlabel ('t[s]'); pause
plot(t, output_signal); title ('Output Signal - sum of modulated I and Q'); ylabel ('Amplitude');  xlabel ('t[s]'); pause;



% transmittiing through AWGN

%calculating SNR
snr = EbNo + 10*log10(k);

noise=randn(size(output_signal)); % random noise generation
constant=std(output_signal)/(std(noise)*10^(snr/20));
received_signal=output_signal + noise*constant; %output of transmitter


%demodulation

I_recovered=received_signal.*carrier_I;
Q_recovered=received_signal.*carrier_Q;


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
mappedSignal=mapped(:,2:3) %just to display and compare with receivedSignal

%plots, constellation
fig = scatterplot(complex_constell); pause;
hold on;
scatterplot(mapped(:,2:3), 1, 0,'ro',fig); pause;
scatterplot(receivedSignal, 1, 0, 'g.', fig); pause;
hold off;
close all;

%inverse mapping and symbols to bitstream conversion
recovered_constellation=rec_constell(receivedSignal, map);
recovered_symbols=gray2symbols(recovered_constellation, map); 
recovered_bits=symbols2bits(recovered_symbols);

%plots in frequency most importatnt blue and red!
hold on;
plot_fft(sI, fs, 'b', 'Spectrum of modulating signal sI'); pause
plot_fft(modulated_I, fs, 'g', 'Spectrum of modulated signal (sI*carrier)'); pause
plot_fft(I_recovered, fs, 'y', 'Spectrum of modulated signal after multiplication by carrier at the receiver'); pause
plot_fft(2*filtered_I, fs, 'r', 'Spectrum of demodulated signal after LPF and amplitude compensation'); pause
hold off;

%plots time
plot(t, sI, 'b'); title('Modulating signal I'); xlabel('t[s]'); ylabel('A'); pause
hold on;
plot(t, filtered_I, 'r'); title('Recovered I'); xlabel('t[s]'); ylabel('A'); pause
hold off;
%plot(t, modulated_I, 'r'); title('Modulated signal (sI*carrier)');xlabel('t[s]'); ylabel('A'); pause
%hold on;
%plot(t, I_recovered, 'g'); title('Modulated signal after multiplication by carrier at the receiver');xlabel('t[s]'); ylabel('A'); pause

plot(t, sQ, 'b'); title('Modulating signal Q'); xlabel('t[s]'); ylabel('A'); pause
hold on;
plot(t, filtered_Q, 'r'); title('Recovered Q'); xlabel('t[s]'); ylabel('A'); pause

%real ber
[bit_errors, ber]=ber_calc(recovered_bits, stream);

%theoretical ber
x=sqrt(3*k*EbNo/(M-1));
theoretical_ber=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));

close all;
bit_errors
ber
theoretical_ber
length(recovered_bits)
