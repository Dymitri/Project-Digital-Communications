EbNodB=-6:2:24
EbNo=10.^(EbNodB/10);
k=8;
M=2^k;
x=sqrt(3*k*EbNo/(M-1));
Pb=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));
semilogy(EbNodB,Pb)