function b = RandBitStream(n)
%generates random bitstream of length n
    x = randint(n,1,[1,0]);
    b = x.';
end 