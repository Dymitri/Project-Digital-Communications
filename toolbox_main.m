M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 30000;                  % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor


rng('default')              % Use default random number generator
dataIn = randi([0 1],n,1);  % Generate vector of binary data
%%dataIn_ok = RandBitStream(n); % THE SAME??!!

% plotting first 40 bits, needed?

%stem(dataIn(1:40),'filled');
%title('Random Bits');
%xlabel('Bit Index'); ylabel('Binary Value');


% changing bits to ints, needed for toolbox qam function

dataInMatrix = reshape(dataIn, length(dataIn)/4, 4); % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers

% ploting ints

figure; % Create new figure window.
stem(dataSymbolsIn(1:10));
title('Random Symbols');
xlabel('Symbol Index'); ylabel('Integer Value');


%modulation with use of one function

dataMod = qammod(dataSymbolsIn, M);


%calculating SNR

EbNo = 10;
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

%transmittiing through AWGN

receivedSignal = awgn(dataMod, snr, 'measured');

%creating constellation diagram

sPlotFig = scatterplot(receivedSignal, 1, 0, 'g.');
hold on
scatterplot(dataMod, 1, 0, 'k*', sPlotFig);

%demodulating using one function

dataSymbolsOut = qamdemod(receivedSignal, M);

%symbol to bit mapping

dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:);                 % Return data in column vector

%BER

[numErrors, ber] = biterr(dataIn, dataOut);
fprintf('\nThe bit error rate = %5.2e, based on %d errors\n', ...
    ber, numErrors)