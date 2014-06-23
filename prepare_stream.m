function [ stream, len ] = prepare_stream( bitstr, k )
%in case of wrong bitstream length
%zero padding

len=length(bitstr);

while mod(len,k) ~= 0
    if k>len
        k
        len
        pads=k-len
    else
        pads=k*ceil(len/k)-len;
    end
    
   stream=padarray(bitstr,[0 pads], 'post');
   len=length(stream)
   pause;
end

