
 clc; clear; close all;

  M = 16;
 n = 1000; % number of bits in bitstream
 k=log2(M);
 M = 16; % order of the modulation
 n = 20; % number of bits in bitstream
 k=log2(M); % size of one symbol?
 EbNo = 10; % to calculate snr
 f=10; %frequency of the carrier
 T=0.1; %symbol time ?
 fsym=1/T; %freq of symbols;% rolloff = 0.5;
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

[I, Q]=split_stream(mapped);




%fig=plot(real(complex_constell),imag(complex_constell),'b.');
t=0:1:length(I)-1;
fig = scatterplot(complex_constell);
hold on;
pause;
scatterplot(mapped(:,2:3), 1, 0,'ro',fig);
pause;

% 
%sI=I/2.*square(2*pi*fsym*t)'+(I/2);
%sQ=Q/2.*square(2*pi*fsym*t)'+(Q/2);
nn=length(I);

t=0:.01:nn;
x=1:1:(nn+1)*100;
for i=1:nn
    for j=1:.1:i+1;
        sI(x(i*100:(i+1)*100))=I(i);
        sQ(x(i*100:(i+1)*100))=Q(i);
    end
end
sI=sI(100:end);
sQ=sQ(100:end);
hold off;
plot(t, sI, 'r');
pause
plot(t, sQ, 'g');
pause

%t=0:step:1;
carrier_I=cos(2*pi*f*t);%Carrier 
carrier_Q=-sin(2*pi*f*t);


modulated_I=sI.*carrier_I;
modulated_Q=sQ.*carrier_Q;

output_signal=modulated_I+modulated_Q;

plot(t, modulated_I);
pause;
plot(t, modulated_Q);
pause
plot(t, output_signal);
pause;


snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol); %%?

% %transmittiing through AWGN
% 
% receivedSignal = awgn(mapped(:,2:3), snr, 'measured'); %% to be replaced
%
% output_signal + noise !!!!!!!!!!!!!!!!
noise=randn(size(mapped(:,2:3))); % random noise generation
constant=std(mapped(:,2:3))/(std(noise)*10^(snr/20));
receivedSignal=mapped(:,2:3) + noise*constant; %output of transmitter
noise1=noise*constant;



%replace outpu_signal to that with noise;
I_recovered=output_signal'*carrier_I;
Q_recovered=output_signal'*carrier_Q;

plot(t, I_recovered);
pause


input=I_recovered;
fs=8*f;
dt=1/fs;
N=length(input);
df = 1/(N*dt);
ff = df * (0 : N/2);
Y=fft(input)/N;
YY=Y(1:((N+1)/2));
plot(ff,abs(YY));

[Y, fff]=fft_ok(I_recovered, 8*f);
plot(fff, abs(Y));





scatterplot(receivedSignal, 1, 0, 'g.', fig);

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
