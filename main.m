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

  
 M = 256; % order of the modulation
 n = 100; % number of bits in bitstream
 k=log2(M); % bits per symbol
 EbNo = 20; % SNR per bit
 fc=4*n; %frequency of the carrier
 snr = EbNo + 10*log10(k); %signal to noise ratio
 fs = 10*fc;          %desired sampling frequency

 
stream=RandBitStream(n); %generating random bitstream
[stream, nb]=prepare_stream(stream, k); %zero padding if needed

frs=nb*round(fs/nb); % real sampling frequency [Hz]
Ts=1/frs;	        % Sampling time [s]
Tb=1/nb;		        % Bit duration time [s]
t=[0:Ts:(nb*Tb)-Ts]';% Time vector initialization


% Modulator

%mapping
[ map, complex_constell ]=gray(M); %creating constellation
grouped_bits=bitgrp(stream, M); %grouping bits
symbols=binary2dec(grouped_bits); % creating symbols
mapped=map2gray(symbols, map); %mapping symbols to the constellation points



%splitting
[I, Q]=split_stream(mapped);



% creating modulating signals sI and sQ

len=nb/k; %length of I or Q in symbols
rep=frs/len; %number of repetitions of a single bit (to obtain square wave)


for i=1:len
        sI(i*rep:(i+1)*rep)=I(i);
        sQ(i*rep:(i+1)*rep)=Q(i);
end
sI=sI(end-frs:end-1); %modulating signals
sQ=sQ(end-frs:end-1);



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




received_signal=add_noise(output_signal, snr);
%received_signal = awgn(output_signal, snr, 'measured');

%demodulation

I_recovered=received_signal.*carrier_I;
Q_recovered=received_signal.*carrier_Q;


%filtration

filtered_I=lpf(I_recovered, fc, frs);
filtered_Q=lpf(Q_recovered, fc, frs);


filtered_I=2*filtered_I;
filtered_Q=2*filtered_Q;

%averaging


demodulated_I=mean(reshape(filtered_I, ceil(frs/len), []))';
demodulated_Q=mean(reshape(filtered_Q, ceil(frs/len), []))';
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
plot_fft(sI, frs, 'b', 'Spectrum of modulating signal sI'); pause
plot_fft(modulated_I, frs, 'g', 'Spectrum of modulated signal (sI*carrier)'); pause
plot_fft(I_recovered, frs, 'y', 'Spectrum of modulated signal after multiplication by carrier at the receiver'); pause
plot_fft(2*filtered_I, frs, 'r', 'Spectrum of demodulated signal after LPF and amplitude compensation'); pause
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

[th_ber,th_ber_haykin, pr_ber, snr_values]= ber_data(output_signal, -20 , 20, 1, M, k, carrier_I, carrier_Q, fc, frs, EbNo, snr, len, map, mapped, stream);
close all;
bit_errors
ber
theoretical_ber



semilogy(snr_values,th_ber,'bo-');
hold on;
semilogy(snr_values,th_ber_haykin,'r.-');
semilogy(snr_values,pr_ber,'mx-');
%axis([-30 10 10^-5 0.5]) 
grid on
legend('theory', 'theory_haykin', 'simulation');
xlabel('SNR, dB');
ylabel('Bit Error Rate');

length(recovered_bits) 
