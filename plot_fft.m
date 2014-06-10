function [ output_args ] = plot_fft( signal, f, c, plot_title )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[Y, fff]=fft_ok(signal, f);
plot(fff, abs(Y), c); title (plot_title); xlabel ('f [Hz]'); 
ylabel ('|A|');

end

