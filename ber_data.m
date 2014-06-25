function [ th_ber, th_ber_haykin, pr_ber, bit_err, snr_vald ] = ber_data( output_signal, start_snr, stop_snr, step, M, k, carrier_I, carrier_Q, frs, h_lpf, len, map, mapped, stream)
%calculates vectors used to compare theoretial bit error level vs practical

snr_vald=start_snr:step:stop_snr;
snr_val=10.^(snr_vald/20);
snr_val=snr_val./0.6;

for i=1:1:length(snr_val)
    
    received_signal=add_noise(output_signal, snr_vald(i));

  

    %demodulation
    I_recovered=received_signal.*carrier_I;
    Q_recovered=received_signal.*carrier_Q;


    %filtration    
    filtered_I=filter(h_lpf,I_recovered);
    filtered_Q=filter(h_lpf,Q_recovered);

    filtered_I=2*filtered_I;
    filtered_Q=2*filtered_Q;

    %averaging
    demodulated_I=mean(reshape(filtered_I, ceil(frs/len), []))';
    demodulated_Q=mean(reshape(filtered_Q, ceil(frs/len), []))';
    receivedSignal=horzcat(demodulated_I, demodulated_Q);
    mappedSignal=mapped(:,2:3); %just to display and compare with receivedSignal


    %inverse mapping and symbols to bitstream conversion
    recovered_constellation=rec_constell(receivedSignal, map);
    recovered_symbols=gray2symbols(recovered_constellation, map); 
    recovered_bits=symbols2bits(recovered_symbols);



    %real ber
    [bit_errors, ber]=ber_calc(recovered_bits, stream);

    bit_err(i)=bit_errors;
    pr_ber(i)=ber;
end



x=sqrt(3*k*snr_val/(M-1));
th_ber=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));
th_ber_haykin=0.5*(1-(1/sqrt(M)))*erfc(sqrt(snr_val));

hold off;


end

