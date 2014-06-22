function [ output_args ] = plot_fft( signal, f, c, plot_title )
%plots spectrum of the given signal
%
%signal - input signal in time domain
%f - sampling frequency
%c - color of the plot
%plot_title - title of the plot

[Y, fff]=fft_ok(signal, f);
plot(fff, abs(Y), c); title (plot_title); xlabel ('f [Hz]'); 
ylabel ('|A|');

end

