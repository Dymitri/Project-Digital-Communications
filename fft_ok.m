function [ YY, ff ] = fft_ok( signal, fs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt=1/fs;
N=length(signal);
df = 1/(N*dt);
ff = df * (0 : N/2-1);
Y=fft(signal)/N;
YY=Y(1:(N/2));


end

