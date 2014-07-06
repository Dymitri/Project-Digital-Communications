function [ th_ber, th_ber_haykin, pr_ber, bit_err, snr_vald ] = ber_data( n, start_snr, stop_snr, step, M, k, carrier_I, carrier_Q, frs, fc, h_lpf, len, map, mapped, stream)
%calculates vectors used to compare theoretial bit error level vs practical

snr_vald=start_snr:step:stop_snr;
snr_val=10.^(snr_vald/20);
snr_val=snr_val./0.6;



stream=RandBitStream(n); %generating random bitstream
[stream, nb]=prepare_stream(stream, k); %zero padding if needed


%frs=nb*ceil(fs/nb); % real sampling frequency [Hz]
Ts=1/frs;	        % Sampling time [s]
Tb=1/nb;		        % Bit duration time [s]
t=[0:Ts:(nb*Tb)-Ts]';% Time vector initialization

% Modulator

%mapping
[ map, complex_constell ]=gray(M); %creating constellation
grouped_bits=bitgrp(stream, M); %grouping bits
symbols=binary2dec(grouped_bits); % creating symbols
mapped=map2gray(symbols, map); %mapping symbols to the constellation points



%splitting
[I, Q]=split_stream(mapped);



% creating modulating signals sI and sQ

len=nb/k; %length of I or Q in symbols
rep=frs/len; %number of repetitions of a single bit (to obtain square wave)


for i=1:len  
        sI(i*rep:(i+1)*rep)=I(i);
        sQ(i*rep:(i+1)*rep)=Q(i);
end
sI=sI(end-frs:end-1); %modulating signals
sQ=sQ(end-frs:end-1);



% Carriers
carrier_I=cos(2*pi*fc*t);
carrier_Q=-sin(2*pi*fc*t);

% Modulating
modulated_I=sI'.*carrier_I;
modulated_Q=sQ'.*carrier_Q;

% Summing I and Q
output_signal=modulated_I+modulated_Q;



% transmittiing through AWGN



for i=1:1:length(snr_val)

received_signal=addnoise(output_signal, 10*snr_vald(i), Tb, rep);
%received_signal = awgn(output_signal, snr, 'measured');

%demodulation

I_recovered=received_signal.*carrier_I;
Q_recovered=received_signal.*carrier_Q;


%filtration
%h_lpf=lpf(fc, frs, filter_order);

filtered_I=filter(h_lpf,I_recovered);
filtered_Q=filter(h_lpf,Q_recovered);




filtered_I=2*filtered_I;
filtered_Q=2*filtered_Q;

%averaging


demodulated_I=mean(reshape(filtered_I, ceil(frs/len), []))';
demodulated_Q=mean(reshape(filtered_Q, ceil(frs/len), []))';
receivedSignal=horzcat(demodulated_I, demodulated_Q)
mappedSignal=mapped(:,2:3) %just to display and compare with receivedSignal


%inverse mapping and symbols to bitstream conversion
recovered_constellation=rec_constell(receivedSignal, map);
recovered_symbols=gray2symbols(recovered_constellation, map); 
recovered_bits=symbols2bits(recovered_symbols);




%real ber
[bit_errors, ber]=ber_calc(recovered_bits, stream);
    pr_ber(i)=ber;
    bit_err(i)=bit_errors;
end



x=sqrt(3*k*snr_val/(M-1));
th_ber=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));

th_ber_haykin=0.5*(1-(1/sqrt(M)))*erfc(sqrt(snr_val));

plot(I_recovered);

end

