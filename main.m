
 clc; clear; close all;

 M = 16;
 n = 100; % number of bits in bitstream
 k=log2(M)
 
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




fig=scatterplot(complex_constell);
hold on;
pause;
scatterplot(mapped(:,2:3), 1, 0,'ro',fig);
pause;


EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol); %%?
%transmittiing through AWGN

receivedSignal = awgn(mapped(:,2:3), snr, 'measured'); %% to be replaced


scatterplot(receivedSignal, 1, 0, 'g.', fig);

recovered_constellation=rec_constell(receivedSignal);

scatterplot((2.*round((recovered_constellation+1)/2)-1), 1, 0, 'ko', fig);




