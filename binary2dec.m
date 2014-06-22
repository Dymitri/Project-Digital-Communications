function [ outvec ] = binary2dec( invec )
%convers vectors of bits to the decimal numbers

outvec = invec*(2.^(size(invec,2)-1:-1:0))';

end

