
 clc; clear; close all;

 M = 16;
 n = 1000; % number of bits in bitstream
 k=log2(M);
 
% rolloff = 0.5;
% snr = 20;
% kanal = 2;
% BR = 64000;

numSamplesPerSymbol = 1;    % Oversampling factor
%Tb = 1/BR;     %perioda 1 bit
%tau = 1.0e-6;
%shift = round(tau/Tb);

stream=RandBitStream(n);


% Modulator

%mapping
[ map, complex_constell ]=gray(M); %constellation


grouped_bits=bitgrp(stream, M);
symbols=binary2dec(grouped_bits);
mapped=map2gray(symbols, map);

%splitting

%[I, Q]=split_stream(mapped);




fig=plot(real(complex_constell),imag(complex_constell),'b.');
hold on;
pause;
scatterplot(mapped(:,2:3), 1, 0,'ro',fig);
pause;

% 
EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol); %%?
% %transmittiing through AWGN
% 
% receivedSignal = awgn(mapped(:,2:3), snr, 'measured'); %% to be replaced
noise=randn(size(mapped(:,2:3))); % random noise generation
constant=std(mapped(:,2:3))/(std(noise)*10^(snr/20));
receivedSignal=mapped(:,2:3) + noise*constant; %output of transmitter
noise1=noise*constant;

scatterplot(receivedSignal, 1, 0, 'g.', fig);

recovered_constellation=rec_constell(receivedSignal, map);

scatterplot((2.*round((recovered_constellation+1)/2)-1), 1, 0, 'ko', fig);

recovered_symbols=gray2symbols(recovered_constellation, map);
recovered_bits=symbols2bits(recovered_symbols);

[bit_errors, ber]=ber_calc(recovered_bits, stream);

bit_errors
ber
length(recovered_bits)
