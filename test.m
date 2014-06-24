SNR = 5; frameLen = 100;
x = rand(100, 1) > 0.5;
y = awgn(2*x-1, SNR);
z = y > 0;
biterr(x, z)