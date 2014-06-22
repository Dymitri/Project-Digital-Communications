function [ YY, ff ] = fft_ok( signal, fs )
%calculates frequency spectrum of the given signal
%
%signal -  input signal in time domain
%fs - sampling frequency
%YY - amplitude of the output signal
%ff - values of the frequency axis

dtt=1/fs;
N=length(signal);
df = 1/(N*dtt);
ff = df * (0 : N/2-1);
Y=fft(signal)/N;
YY=Y(1:(N/2));


end

