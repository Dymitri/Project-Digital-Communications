function [ output_args ] = plot_fft( signal, f, c )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Y, fff]=fft_ok(signal, 8*f);
plot(fff, abs(Y), c);

end

