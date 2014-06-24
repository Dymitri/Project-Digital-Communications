function [ th_ber, th_ber_haykin, pr_ber, snr_val ] = ber_data( output_signal, start_snr, stop_snr, step, M, k, carrier_I, carrier_Q, fc, frs, EbNo, snr, len, map, mapped, stream)
%calculates vectors used to compare theoretial bit error level vs practical

snr_val=start_snr:step:stop_snr;
snr_val=10.^(snr_val/10);

for i=1:1:length(snr_val)
    
    snr_val=10.^(snr_val/10);
    received_signal=add_noise(output_signal, snr_val(i));

  

    %demodulation
    I_recovered=received_signal.*carrier_I;
    Q_recovered=received_signal.*carrier_Q;


    %filtration
    filtered_I=lpf(I_recovered, fc, frs);
    filtered_Q=lpf(Q_recovered, fc, frs);


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

    %theoretical ber
    %x=sqrt(3*k*EbNo/(M-1));
    x=sqrt(3*k*snr_val(i)/(M-1));
    theoretical_ber=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));
    theoretical_ber_haykin=0.5*(1-(1/sqrt(M)))*erfc(sqrt(snr_val(i)));
    
    pr_ber(i)=ber;
    th_ber(i)=theoretical_ber;
    th_ber_haykin(i)=theoretical_ber_haykin;
end

end

