function [ n, ber ] = ber_calc( received, sent )
%compares generated bitstream and received one
%
%n - number of flipped bits
%ber - bit error rate
n=0;
for i=1:length(received)
    if received(i)~=sent(i)
    n=n+1;    
    end
    
end
ber=n/length(received);
end
