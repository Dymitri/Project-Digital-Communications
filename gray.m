function r = gray(num)
    


o=bitxor(bitshift(num,1),num);

% return (num >> 1) ^ num;
i=0;
while 2^i < num
    i=i+1;
end

r=bitand(o, 3);
i
num
end 