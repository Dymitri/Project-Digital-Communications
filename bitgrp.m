function symb = bitgrp(bitstr, m)
len=length(bitstr);
bitspersymbol=log2(m);
i=1;
j=1;
%while i < (len/(bitspersymbol+1))
while i< 100
for k = (i*bitspersymbol) : (2*i*bitspersymbol)
   symb(j, i)=bitstr(k, : )*2^(i-1);
    i=i+1;
    end
    k=k+1;
end
k
j=j+1;
end

%Y = typecast(X, type)



function symb = bitgrp(bitstr, m)
len=length(bitstr)
bps=log2(m)

while i < 