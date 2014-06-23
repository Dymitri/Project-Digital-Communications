function [ stream, len ] = prepare_stream( bitstr, k )
%in case of wrong bitstream length
%zero padding
len=length(bitstr);

    if mod(len,k) == 0
        pads=0;
    end


    if mod(len,k) ~= 0
        if k>len
            pads=k-len;
        else
            pads=k*ceil(len/k)-len;
        end
    end
    
        
    
        stream=padarray(bitstr,[0 pads], 'post');
        len=length(stream);
    end