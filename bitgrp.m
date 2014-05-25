function symb = bitgrp(bitstr, m)
len=length(bitstr);
bitspersymbol=log2(m);
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
