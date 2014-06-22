function [map mod] = gray(M)
% returns map of gray coded points in a form: decimal value, I, Q
% and mod in a complex form: I+jQ

k = log2(M); % number of bits in each constellation


alphaRe = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
alphaIm = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
% input - decimal equivalent of all combinations with b0b1b2b3

ip = [0:(M-1)];
ipBin = dec2bin(ip.'); % decimal to binary
% taking b0b1 for real
ipDecRe = bin2dec(ipBin(:,[1:k/2]));
ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
% taking b2b3 for imaginary
ipDecIm = bin2dec(ipBin(:,[k/2+1:k]));
ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2));
% mapping the Gray coded symbols into constellation
modRe = alphaRe(ipGrayDecRe+1);
modIm = alphaIm(ipGrayDecIm+1);
% complex constellation
mod = modRe + j*modIm;

map=[ip; modRe; modIm]';
end