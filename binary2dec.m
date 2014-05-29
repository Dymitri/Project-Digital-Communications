function [ outvec ] = binary2dec( invec )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
outvec = invec*(2.^(size(invec,2)-1:-1:0))';

end

