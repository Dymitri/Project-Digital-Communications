
 clc; clear; close all;

 M = 16;
 n = 100; % number of bits in bitstream
 
% rolloff = 0.5;
% snr = 20;
% kanal = 2;
% BR = 64000;

%nsamp = 4; % Oversampling rate
%Tb = 1/BR;     %perioda 1 bit
%tau = 1.0e-6;
%shift = round(tau/Tb);

stream=RandBitStream(n);


% Modulator

%mapping
[ map, complex_constell ]=gray(M); %constellation
scatterplot(complex_constell);

grouped_bits=bitgrp(stream, M);
symbols=binary2dec(grouped_bits);
mapped=map2gray(symbols, map);

%splitting

[I, Q]=split_stream(mapped);



%now we have to map symbols to points in map.
%map2gray function is not finished yet

