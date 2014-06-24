function [ nsignal ] = add_noise( signal, snr )
% awgn


noise=randn(size(signal)); % random noise generation
constant=std(signal)/(std(noise)*10^(snr/20));
nsignal=signal + noise*constant; %input of the receiver



end

