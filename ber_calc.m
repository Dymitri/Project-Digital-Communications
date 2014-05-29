function [ n, ber ] = ber_calc( received, sent )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n=0;
for i=1:length(received)
    if received(i)~=sent(i)
    n=n+1;    
    end
    
end
ber=n/length(received);
end
