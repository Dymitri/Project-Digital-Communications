
 clc; clear; close all;

 M = 16;
 n = 100; % number of bits in bitstream
 
% rolloff = 0.5;
% snr = 20;
% kanal = 2;
% BR = 64000;

nsamp = 4; % Oversampling rate
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

[I, Q]=split_stream(mapped);



% Filter Definition
% Define filter-related parameters.
filtorder = 40; % Filter order
delay = filtorder/(nsamp*2); % Group delay (# of input samples)
 rolloff = 0.25; % Rolloff factor of filter

% Create a square root raised cosine filter.
rrcfilter = rcosine(1,nsamp,'fir/sqrt',rolloff,delay);

% Transmitted Signal
% Upsample and apply square root raised cosine filter.
ytx = rcosflt(mapped,1,nsamp,'filter',rrcfilter);
fig=scatterplot(complex_constell);
hold on;
pause;
scatterplot(ytx(:,2:3), 1, 0,'g.',fig);




%now we have to map symbols to points in map.
%map2gray function is not finished yet

