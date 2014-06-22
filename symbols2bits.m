function [ rec_bitstream ] = symbols2bits( symbols )
% symbols to bits conversion


matrix=dec2bin(symbols);
str=reshape(matrix', [], 1)';

for i= 1:length(str)
    
    rec_bitstream(i)=str2num(str(i));
    
end

