function symb = bitgrp(bitstr, m)
% returns array of a grouped bits of lengrh depending on the modulation
% order m used

len=length(bitstr);
bitspersymbol=log2(m)

%in case of wrong bitstream length
%zero padding
num_sym=len/bitspersymbol;
while mod(len,bitspersymbol) ~= 0
    pads=bitspersymbol*ceil(num_sym)- len
    if len<=bitspersymbol
        pads=pads+bitspersymbol;
    end
    bitstr=padarray(bitstr,[0 pads], 'post')
    len=length(bitstr)
    num_sym=len/bitspersymbol
    modd=mod(num_sym,bitspersymbol) 
    pause
end
    
    
i=1;
j=1;
k=1;
while k<(len+1)
    while i<(bitspersymbol+1)
        symb(j,i)=bitstr(k);
        k=k+1;
        i=i+1;
    end
i=1;
j=j+1;
end
