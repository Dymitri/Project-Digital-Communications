function [ Hd ] = lpf( fp, fs, n_lpf )
%n_lpf order Butterworth low pass filter
%
%fp - cutoff frequency
%fs - sampling frequency
d=fdesign.lowpass('N,Fc',10,fp,fs);
designmethods(d);

Hd = design(d,'butter');

clc;
end

