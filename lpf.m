function [ filtered_signal ] = lpf( signal, fp, fs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%d=fdesign.lowpass('Fp,Fst,Ap,Ast',0.15,0.25,1,60);

%Create a filter of order 10 with a 6-dB frequency of 9.6 kHz and a sampling frequency of 48 kHz.

d=fdesign.lowpass('N,Fc',10,fp,fs);
designmethods(d);

Hd = design(d,'butter');
%fvtool(Hd);

filtered_signal = filter(Hd,signal);
clc;
end

