function [ filtered_signal ] = lpf( signal, fp, fs )
%10 order Butterworth low pass filter
%
%signal - input signal in time domain
%fp - cutoff frequency
%fs - sampling frequency
%filtered_signal - output signal in time domain
d=fdesign.lowpass('N,Fc',10,fp,fs);
designmethods(d);

Hd = design(d,'butter');

filtered_signal = filter(Hd,signal);
clc;
end

