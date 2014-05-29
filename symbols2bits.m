function [ rec_bitstream ] = symbols2bits( symbols )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


matrix=dec2bin(symbols);
str=reshape(matrix', [], 1)';
%rec_bitstream=str2num(str);

for i= 1:length(str)
    
    rec_bitstream(i)=str2num(str(i));
    
end

